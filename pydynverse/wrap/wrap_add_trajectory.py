import pandas as pd
import pandas.api.types as pdtypes

from .._logging import logger
from .convert_milestone_percentages_to_progressions import convert_milestone_percentages_to_progressions
from .convert_progressions_to_milestone_percentages import convert_progressions_to_milestone_percentages
from .classify_milestone_network import classify_milestone_network

def add_trajectory(
    dataset:dict,
    milestone_network:pd.DataFrame,
    milestone_ids:list=None,
    divergence_regions:pd.DataFrame=None,
    milestone_percentages:pd.DataFrame=None,
    progressions:pd.DataFrame=None,
    allow_self_loops=False,
    **kwargs,
):
    cell_ids = dataset["cell_ids"]

    if milestone_ids is None:
        # 如果milestone_ids没有给定的话, 从milestone_network的from和to中来唯一获得
        milestone_ids = list(set(milestone_network["from"]) | set(milestone_network["to"]))

    # 里程碑网络结构
    milestone_network = check_milestone_network(milestone_ids, milestone_network, allow_self_loops)

    # 延迟承诺区域
    if divergence_regions is None:
        divergence_regions = None
    divergence_regions = check_divergence_regions(milestone_ids, divergence_regions)

    # 检测milestone percentages和progressions
    if  (milestone_percentages is None) == (progressions is None):
        # 同时为空或同时不为空
        if not milestone_percentages is None:
            logger.warning("Both milestone_percentages and progressions are given, will only use progressions")
            milestone_percentages = None
        else:
            error_msg = "Exactly one of milestone_percentages or progressions, must be defined, the other must be None"
            logger.error(error_msg)
            raise ValueError(error_msg)
       

    if progressions is None:
        # milestone_percentages -> progressions, addtrajectory测试用例执行这里
        check_milestone_percentages(cell_ids, milestone_ids, milestone_percentages)
        progressions = convert_milestone_percentages_to_progressions(cell_ids, milestone_ids, milestone_network, milestone_percentages)
        progressions = check_progressions(cell_ids, milestone_ids, milestone_network, progressions)
    else:
        # TODO: progressions -> milestone_percentages, slingshot执行这里, add_branch_trajectory测试用例执行这里
        progressions = check_progressions(cell_ids, milestone_ids,milestone_network, progressions)
        milestone_percentages = convert_progressions_to_milestone_percentages(cell_ids, milestone_ids, milestone_network, progressions)
        milestone_percentages = check_milestone_percentages(cell_ids, milestone_ids, milestone_percentages)

    # TODO: 检测延迟承诺区域是否把对应的细胞覆盖了

    # 添加输出结构
    trajectory_dict = {
        "pydynwrap::with_trajectory": True,
        "milestone_ids": milestone_ids,
        "milestone_network": milestone_network,
        "divergence_regions": divergence_regions,
        "milestone_percentages": milestone_percentages,
        "progressions": progressions,
    }
    dataset.update(trajectory_dict)
    dataset["pydynwrap:with_trajectory"] = True

    # 检测milestone网络拓扑的结构类型
    dataset = changed_topology(dataset)

    return dataset


def is_wrapper_with_trajectory(trajectory):
    return trajectory.get("pydynwrap:with_trajectory", False)



def check_milestone_network(milestone_ids, milestone_network, allow_self_loops):
    # 判断列名合法性
    assert set(["from", "to", "length", "directed"]).issubset(set(milestone_network.columns))
    
    milestone_network = milestone_network[["from", "to", "length", "directed"]]

    # 判断各列的元素类型
    assert pdtypes.is_string_dtype(milestone_network["from"])
    assert pdtypes.is_string_dtype(milestone_network["to"])
    assert pdtypes.is_numeric_dtype(milestone_network["length"])
    assert pdtypes.is_bool_dtype(milestone_network["directed"])
    # 判断from和to两列的milestone_id是否合法
    assert set(milestone_network["from"]).issubset(set(milestone_ids))
    assert set(milestone_network["to"]).issubset(set(milestone_ids))
    # 判断from和to两列的行元素不能重复
    assert (not milestone_network[["from", "to"]].duplicated().any())

    if not allow_self_loops:
        # 判断不允许自环时，每行的from和to两列的值不能重复
        assert (not (milestone_network["from"]==milestone_network["to"]).any())

    check = milestone_network
    if allow_self_loops:
        # 过滤掉自环(一跳环)的milestone
        check = check[check["from"]==check["to"]]
    check = pd.merge(check, check, left_on="to", right_on="from")
    # 判断二跳环
    assert (not (check["from_x"]==check["to_y"]).any())

    return milestone_network


# 这里的一堆assert暂时不用管
def check_divergence_regions(milestone_ids, divergence_regions):
    return divergence_regions


def check_milestone_percentages(cell_ids, milestone_ids, milestone_percentages):
    return milestone_percentages


def check_progressions(cell_ids, milestone_ids, milestone_network, progressions):
    return progressions


def changed_topology(trajectory):
    classification = classify_milestone_network(trajectory["milestone_network"])
    trajectory["trajectory_type"] = classification["network_type"]
    trajectory["directed"] = classification["directed"]
    return trajectory
