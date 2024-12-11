import numpy as np
import pandas as pd

from .wrap_add_expression import get_expression


def add_dimred(
        dataset,
        dimred,
        dimred_milestones=None,  # 里程碑
        dimred_segment_progressions=None,  # 不连通的分段过程
        dimred_segment_points=None,  # 不连通的分段点
        project_trajectory=True,  # 投影轨迹
        connect_segments=False,  # 是否对分段进行连接
        expression_source="expression"):
    # 添加降维结果

    cell_ids = dataset["cell_ids"]
    dimred_type = type(dimred)
    if dimred_type == np.array:
        # 指定降维结果,直接赋值
        dimred = process_dimred(dataset, dimred)
    else:
        # 指定降维方法, 需要传入表达矩阵执行来获得降维结果
        expression = get_expression(dataset, expression_source)
        dimred = dimred(expression)
        dimred.index = cell_ids  # 这里的DataFrame手动添加索引

    if dimred_milestones:
        # TODO: 降维后的里程碑milestone设置
        pass

    if dimred_segment_points or dimred_segment_progressions:
        # TODO: 降维后的不连通的分段设置
        pass

    dataset["dimred"] = dimred
    dataset["dimred_milestones"] = dimred_milestones
    dataset["dimred_segment_progressions"] = dimred_segment_progressions
    dataset["dimred_segment_points"] = dimred_segment_points

    return dataset


def is_wrapper_with_dimred(dataset):
    if "dimred" in dataset:
        return True
    else:
        return False


def get_dimred(dataset,
               dimred,
               expression_source="expression",
               return_other_dimreds=False):
    # 获得降维结果
    # 执行quickstart里直接来到这里，跳过了R代码中的一堆判断
    if type(dimred) is pd.DataFrame:
        # 传入的dimred为DataFrame
        dimred = dimred
        extra_out = {}
        
    elif is_wrapper_with_dimred(dataset):
        dimred = dataset["dimred"]
        extra_out = {}
        if return_other_dimreds:
            # 里程碑
            if dataset["dimred_milestones"]:
                extra_out["dimred_milestones"] = dataset["dimred_milestones"]
            # 分段过程
            if (not (dataset["dimred_segment_progressions"] is None)) and (not (dataset["dimred_segment_points"] is None)):
                extra_out["dimred_segment_progressions"] = dataset["dimred_segment_progressions"]
                extra_out["dimred_segment_points"] = dataset["dimred_segment_points"]
    
    dimred.columns = [f"comp_{i+1}" for i in range(dimred.shape[1])] # 这里的降维维度序号从1开始

    dimred = dimred.loc[dataset["cell_ids"]]  # 这步对齐是没有必要的

    return dimred, extra_out


def process_dimred(dataset, dimred, identifier="cell_id", has_rownames=True):
    # TODO: 降维结果dimred为Array或DataFrame时，添加行列名，结果都成了DataFrame了
    pd.DataFrame(dimred)
    return dimred
