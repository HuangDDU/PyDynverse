import numpy as np
import pandas as pd
from tqdm import tqdm
from pydynverse.util.sparse_matrix import is_sparse
from pydynverse.feature.fi_methods import fi_ranger_rf_lite
#TODO:(初步)verison0.0.1

def calculate_feature_importances(X, Y, fi_method=None, verbose=False):
    """
    Calculate feature importance scores.
    
    Parameters:
      X: pd.DataFrame or np.ndarray
         包含各特征（列）的数据框。
      Y: pd.DataFrame, np.ndarray or list
         预测变量的数据框或矩阵，行数需与 X 相同。
      fi_method: dict, optional
         特征重要性方法，默认使用 fi_ranger_rf_lite()。
         该方法应返回一个字典，且其中 key 'fun' 对应一个计算特征重要性的函数。
      verbose: bool, optional
         是否输出额外信息。
    
    Returns:
      pd.DataFrame: 包含三列 predictor_id、feature_id 以及 importance 的数据框，
      其中 predictor_id 来自 Y 的列名，feature_id 来自 X 的列名。
    
    Examples:
      >>> X = pd.DataFrame(np.random.rand(25, 10), columns=[f'feature{i}' for i in range(10)])
      >>> Y = pd.DataFrame(np.random.rand(25, 2), columns=[f'pred{i}' for i in range(2)])
      >>> imp_df = calculate_feature_importances(X, Y)
      >>> print(imp_df)
    """
    # 如果没有指定特征重要性方法，则默认使用 fi_ranger_rf_lite()
    if fi_method is None:
        fi_method = fi_ranger_rf_lite()
    
    # 如果 Y 不是 DataFrame，则转换
    if not isinstance(Y, pd.DataFrame):
        Y = pd.DataFrame(Y)
    
    # 如果 Y 或 X 为稀疏矩阵，则转换为常规矩阵（X通常是expression矩阵）
    if is_sparse(Y):
        Y = pd.DataFrame(Y.toarray(), columns=[str(i) for i in range(Y.shape[1])])
    if is_sparse(X):
        X = pd.DataFrame(X.toarray(), columns=[str(i) for i in range(X.shape[1])])
    
    # 确保 X 为 DataFrame
    if not isinstance(X, pd.DataFrame):
        X = pd.DataFrame(X)
    
    #分别计算每个预测因子的重要性得分
    importances_list = []
    # 针对 Y 中每一个预测变量分别计算特征重要性
    for i in tqdm(range(Y.shape[1]), desc='Generating forests', disable=not verbose):
        if verbose:
            print(f"Generating forest {i+1}/{Y.shape[1]}")
        
        # 提取第 i 列作为目标变量 y
        y = Y.iloc[:, i]
        
        # 如果 y 的数据类型为字符串或布尔，则转换为 category 类型
        if y.dtype == object or y.dtype == bool:
            y = y.astype('category')
        
        # 如果 y 中所有值均相同，则直接返回 0 重要性
        if y.nunique() == 1:
            importance = pd.Series(0, index=X.columns)
        else:
            # 调用 fi_method 返回的 'fun' 函数计算特征重要性
            # 要求 fi_method 的计算函数签名为: fun(X, y, verbose=verbose)
            importance_dict = fi_method['fun'](X, y, verbose=verbose)
            # 将返回结果转换为 Series（保证索引与 X 的列名一致）
            importance = pd.Series(importance_dict, index=X.columns)
        
        # 构造当前预测变量的特征重要性 DataFrame
        predictor_name = Y.columns[i] if Y.columns is not None else str(i)
        temp_df = pd.DataFrame({
            'predictor_id': predictor_name,
            'feature_id': importance.index,
            'importance': importance.values
        })
        importances_list.append(temp_df)
    
    # 合并所有预测变量的结果，并按 importance 降序排列
    final_importances = pd.concat(importances_list, ignore_index=True)
    final_importances = final_importances.sort_values(by='importance', ascending=False).reset_index(drop=True)
    
    return final_importances