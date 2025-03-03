import pydynverse as pdv
import pandas as pd
import numpy as np
from sklearn.metrics import jaccard_score
from pandas import DataFrame, Series


def calculate_mapping(dataset, prediction, grouping='milestones', simplify=False):
    """_summary_

    Args:
        dataset (dict): 原始trajectory数据（需要包含细胞信息和里程碑网络信息）
        prediction (dict): 预测trajectory数据（需要包含细胞信息和里程碑网络信息）
        grouping (str, optional): 选择映射方式，‘milestones’按里程碑映射，’branches‘按照分支映射. Defaults to 'milestones'.
        simplify (bool, optional): 是否简化数据集轨迹. Defaults to True.
    """
    # grouping输入错误则抛出异常方便调试
    if grouping not in ['branches', 'milestones']:
        raise ValueError("grouping must be either 'branches' or 'milestones'")

    if dataset is None or prediction is None:
        return {'recovery': 0, 'relevance': 0, 'F1': 0}

    if simplify:
        dataset = pdv.wrap.simplify_trajectory(dataset)
        prediction = pdv.wrap.simplify_trajectory(prediction)

    # 如果映射到分支上，则传入：
    if grouping == 'branches':
        groups_dataset = pdv.wrap.group_onto_trajectory_edges(dataset)
        groups_prediction = pdv.wrap.group_onto_trajectory_edges(prediction)
    # 如果映射到里程碑上，则传入：
    elif grouping == 'milestones':
        groups_dataset = pdv.wrap.group_onto_nearest_milestones(dataset)
        groups_prediction = pdv.wrap.group_onto_nearest_milestones(prediction)

    # Create a DataFrame to compare the two groups(这里的Dataframe必须同时包含聚类集和cell_ids的信息)
    groups_ref = pd.DataFrame(
        {"ref_map": groups_dataset, "ref_cell": dataset["cell_ids"]})
    groups_pre = pd.DataFrame(
        {"pre_map": groups_prediction, "pre_cell": prediction["cell_ids"]})
    # 处理聚类簇(index是聚类簇，后跟簇内细胞集合)
    groups_ref = groups_ref.groupby('ref_map')['ref_cell'].apply(set)
    groups_pre = groups_pre.groupby('pre_map')['pre_cell'].apply(set)
    # 记录每个group聚类簇的数量(似乎暂时用不上)
    cluster_num_for_ref = len(groups_ref.index)
    cluster_num_for_pre = len(groups_pre.index)
    # 计算得到jaccard指数矩阵(DataFrame)
    jaccard_matrix = pd.DataFrame(
        index=groups_ref.index, columns=groups_pre.index)

    # 计算每个簇之间的 Jaccard 相似度
    for ref_map in groups_ref.index:
        for pre_map in groups_pre.index:
            # 计算交集和并集的大小
            intersection_size = len(
                groups_ref[ref_map].intersection(groups_pre[pre_map]))
            union_size = len(groups_ref[ref_map].union(groups_pre[pre_map]))
            # 计算 Jaccard 相似度
            jaccard_similarity = intersection_size / union_size if union_size != 0 else 0
            # 填入 Jaccard 矩阵
            jaccard_matrix.loc[ref_map, pre_map] = jaccard_similarity

    # 得到jaccard_matrix之后，需要计算recovery 和 relevance
    recovery = jaccard_matrix.max(axis=1).mean()
    relevance = jaccard_matrix.max(axis=0).mean()
    # Calculate F1 score
    F1 = 2 * (recovery * relevance) / (recovery +
                                       relevance) if (recovery + relevance) > 0 else 0

    return {'recovery': recovery, 'relevance': relevance, 'F1': F1}


def calculate_mapping_milestones(dataset, prediction, simplify=False):
    # TODO: 基于pdv.wrap.group_onto_nearest_milestones结果计算
    mapping = calculate_mapping(dataset, prediction, 'milestones', simplify)
    return {f"{key}_milestones": value for key, value in mapping.items()}


def calculate_mapping_branches(dataset, prediction, simplify=False):
    # TODO: 基于pdv.wrap.group_onto_trajectory_edges结果计算
    mapping = calculate_mapping(dataset, prediction, 'branches', simplify)
    return {f"{key}_branches": value for key, value in mapping.items()}
