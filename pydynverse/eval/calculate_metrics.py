import networkx as nx
# from .__init__ import metrics # 交叉重复导入会报错
from ..wrap import simplify_trajectory

from .metric_isomorphic import calc_isomorphic
from .metric_flip import calculate_edge_flip
from .metric_him import calculate_him
from .metric_mapping import calculate_mapping_branches, calculate_mapping_milestones


def calculate_metrics(
    dataset,
    model,
    # metrics=metrics["metric_id"],
    metrics=["isomorphic"],
):
    # 一堆指标检查先不管
    summary_dict = {}
    # 简化轨迹
    dataset = simplify_trajectory(dataset)
    if model is not None:
        model = simplify_trajectory(model)

    # TODO: 其他指标
    # milestone相关指标
    net1 = model["milestone_network"]
    net2 = dataset["milestone_network"]
    if "isomorphic" in metrics:
        # 这里与dynverse不同，也用函数实现
        summary_dict["isomorphic"] = calc_isomorphic(net1, net2)
    if "edge_flip" in metrics:
        summary_dict["edge_flip"] = calculate_edge_flip(net1, net2)
    if "him" in metrics:
        summary_dict["him"] = calculate_him(net1, net2)

    # 传参调整
    if "F1_branch" in metrics:
        summary_dict["F1_branch"] = calculate_mapping_branches()

    if "F1_milestone" in metrics:
        summary_dict["F1_milestone"] = calculate_mapping_milestones()

    # 其他指标
    return summary_dict
