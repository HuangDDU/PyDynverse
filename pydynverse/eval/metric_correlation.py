from pydynverse.wrap.calculate_geodesic_distances import calculate_geodesic_distances
import time
import numpy as np
from scipy.stats import spearmanr
import sys
from typing import Dict, Any
from pydynverse.wrap import select_waypoints

#from ..wrap import is_wrapper_with_waypoint_cells#is_wrapper_with_waypoint_cells接口在哪？？

def calc_correlation(dataset, prediction):
    """
    计算数据集与预测模型之间的地理距离相关性
    param dataset: 包含细胞轨迹信息的字典，需有waypoint_cells、cell_ids等字段
    param prediction: 预测模型的输出结果，结构同dataset
    return: 包含相关性指标和时间消耗的字典
    """
    metrics = {}
    
    # 1. 验证数据结构
    """
    if not is_wrapper_with_waypoint_cells(dataset):
        raise ValueError("Dataset must contain waypoint cells")
    if prediction is not None and not is_wrapper_with_waypoint_cells(prediction):
        raise ValueError("Prediction model must contain waypoint cells")
    """
    # 如果没有预测模型，直接返回0相关性
    if prediction is None:
        return {'correlation': 0.0}
    
    # 2. 确保所有预测细胞都在原始数据中
    pred_cell_ids = prediction.get('cell_ids', [])
    if not all(cell in dataset['cell_ids'] for cell in pred_cell_ids):
        missing = set(pred_cell_ids) - set(dataset['cell_ids'])
        raise ValueError(f"Prediction contains unknown cells: {missing}")
    
    # 统一细胞ID顺序（关键步骤）
    dataset['cell_ids'] = sorted(dataset['cell_ids'])
    prediction['cell_ids'] = dataset['cell_ids']  # 强制使用原始数据顺序

    # 3. 合并waypoint细胞
    waypoints = list(set(dataset.get('waypoint_cells', []) + prediction.get('waypoint_cells', [])))
    waypoints=None #这部分先置空，方便测试
    #这里的calculate_geodesic_distances可能的问题？
    #如果waypoints=[]，那么calculate_geodesic_distances会得到一个空的DataFrame

    # 4. 计算地理距离矩阵
    # 数据集部分
    start_time = time.time()
    #这里传入的参数是waypoint_cells
    dataset_dist = calculate_geodesic_distances(dataset, waypoints)
    metrics['time_waypoint_geodesic'] = time.time() - start_time
    
    # 预测模型部分
    start_time = time.time()
    #这里传入的参数是waypoint_cells
    pred_dist = calculate_geodesic_distances(prediction, waypoints)
    metrics['time_pred_geodesic'] = time.time() - start_time
    
    # 5. 处理无限值（替换为最大浮点数）
    # 处理无限值（关键修复）
    max_float = sys.float_info.max
    def replace_inf(df, max_float):
        return df.replace([np.inf, -np.inf], max_float).to_numpy(dtype=np.float64)
    
    dataset_dist = replace_inf(dataset_dist, max_float)
    pred_dist = replace_inf(pred_dist, max_float)
    
    # 6. 验证矩阵维度一致性
    if dataset_dist.shape != pred_dist.shape:
        raise RuntimeError(f"Distance matrix shape mismatch: {dataset_dist.shape} vs {pred_dist.shape}")
    
    # 7. 计算Spearman相关系数
    start_time = time.time()
    
    # 检查数据是否全相同（避免除零）
    if np.unique(dataset_dist).size == 1 or np.unique(pred_dist).size == 1:
        corr = 0.0
    else:
        # 展平矩阵并计算
        corr, _ = spearmanr(dataset_dist.flatten(), pred_dist.flatten())
        corr = max(corr, 0.0)  # 确保非负
    
    metrics['correlation'] = corr
    metrics['time_correlation'] = time.time() - start_time
    
    return metrics
     