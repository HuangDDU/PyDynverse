import networkx as nx
# from .__init__ import metrics # 交叉重复导入会报错
from ..wrap import simplify_trajectory

from .metric_correlation import calc_correlation
from .metric_isomorphic import calc_isomorphic
from .metric_flip import calculate_edge_flip
from .metric_him import calculate_him
from .metric_mapping import calculate_mapping_branches, calculate_mapping_milestones
import time

def calculate_metrics(
    dataset,
    model,
    simplify=True,
    # metrics=metrics["metric_id"],
    metrics=["isomorphic", "edge_flip"],
):
    # 一堆指标检查先不管
    summary_dict = {}
    # 简化轨迹
    if simplify:
        dataset = simplify_trajectory(dataset)
        model = simplify_trajectory(model)


    if "correlation" in metrics:
        summary_dict["correlation"] = calc_correlation(dataset, model)

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
    """
      # 检查并计算与特征重要性（相关性）相关的指标
    if any(metric in metrics for metric in ["featureimp_cor", "featureimp_wcor"]):
        time0 = time.time()
        featureimp = calculate_featureimp_cor(dataset, model, expression_source=expression_source)
        time1 = time.time()
        summary_dict["time_featureimp"] = time1 - time0
        summary_dict["featureimp_cor"] = featureimp["featureimp_cor"]
        summary_dict["featureimp_wcor"] = featureimp["featureimp_wcor"]

    # 检查并计算与特征重要性（富集分析）相关的指标
    if any(metric in metrics for metric in ["featureimp_ks", "featureimp_wilcox"]):
        time0 = time.time()
        featureimp = calculate_featureimp_enrichment(dataset, model, expression_source=expression_source)
        time1 = time.time()
        summary_dict["time_featureimp_enrichment"] = time1 - time0
        summary_dict["featureimp_ks"] = featureimp["featureimp_ks"]
        summary_dict["featureimp_wilcox"] = featureimp["featureimp_wilcox"]
        """

    # 传参调整
    if "F1_branch" in metrics:
        summary_dict["F1_branch"] = calculate_mapping_branches()

    if "F1_milestone" in metrics:
        summary_dict["F1_milestone"] = calculate_mapping_milestones()

    # 其他指标
    return summary_dict
