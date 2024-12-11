import math
import numpy as np
from .convert_progressions_to_milestone_percentages import convert_progressions_to_milestone_percentages
from .calculate_geodesic_distances import calculate_geodesic_distances


def select_waypoints(
    trajectory,
    n_waypoints=200,
    trafo=math.sqrt,  # R语言中是sqrt函数, 测试的时候写成lambda x:x 更直观
    resolution=None,
    recompute=False
):

    if resolution is None:
        # 如果没有指定resolution分辨率步长，则根据milestone_network网络的总长度与预期waypoint个数来计算
        resolution = trajectory["milestone_network"]["length"].apply(lambda x: trafo(x)).sum()/n_waypoints
    if (not recompute) and (not trajectory.get("waypoints", None) is None):
        # 已经计算过了, 不重复计算
        return trajectory["waypoints"]

    waypoint_progressions = trajectory["milestone_network"].copy()
    waypoint_progressions["percentage"] = waypoint_progressions["length"].apply(lambda x: [i/x for i in np.arange(0, x, resolution)] + [1])  # 按照分辨率为固定步长划分milestone_network的边, 添加到percentage列上
    waypoint_progressions = waypoint_progressions[["from", "to", "percentage"]]
    waypoint_progressions = waypoint_progressions.explode("percentage").reset_index(drop=True)  # 展开
    waypoint_progressions["percentage"] = waypoint_progressions["percentage"].astype("float")

    def waypoint_id_from_progressions_row(row):
        # 从progressions_row的一行综合确定waypoint_id
        match row["percentage"]:
            case 0:
                return f"MILESTONE_BEGIN_W{row['from']}_{row['to']}"
            case 1:
                return f"MILESTONE_END_W{row['from']}_{row['to']}"
            case _:
                return f"W{row.name+1}"  # 这里与R保持一致，W从1开始算
    waypoint_progressions["waypoint_id"] = waypoint_progressions.apply(waypoint_id_from_progressions_row, axis=1)

    waypoint_progressions_tmp = waypoint_progressions.copy()
    waypoint_progressions_tmp = waypoint_progressions_tmp.rename(columns={"waypoint_id": "cell_id"})  # 复用之前列改名
    # 复用milstone的progressions->percentages代码
    waypoint_milestone_percentages = convert_progressions_to_milestone_percentages(
        cell_ids="this argument is unnecessary, I can put everything I want in here!",
        milestone_ids=trajectory["milestone_ids"],
        milestone_network=trajectory["milestone_network"],
        progressions=waypoint_progressions_tmp
    ).rename(columns={"cell_id": "waypoint_id"})

    # waypoint到细胞的距离, (n_waypoints, n_cells)
    waypoint_geodesic_distances = calculate_geodesic_distances(
        trajectory,
        waypoint_milestone_percentages=waypoint_milestone_percentages
    ).loc[waypoint_progressions["waypoint_id"]]

    waypoint_network = waypoint_progressions\
        .sort_values(by=["from", "to", "percentage"])\
        .groupby(["from", "to"])\
        .apply(lambda group: group.assign(
            from_waypoint=group["waypoint_id"],
            to_waypoint=group["waypoint_id"].shift(-1),
        ))\
        .dropna()\
        .reset_index(drop=True)  # 分组后, 组内按照percentage排序, lead函数向当前看下一行的元素, 如果是组内最后一个元素则获得NULL
    waypoint_network = waypoint_network[["from_waypoint", "to_waypoint", "from", "to"]]
    waypoint_network.columns = ["from", "to", "from_milestone_id", "to_milestone_id"]

    waypoints = waypoint_milestone_percentages.iloc[waypoint_milestone_percentages.groupby("waypoint_id")["percentage"].idxmax()].reset_index(drop=True)
    waypoints["milestone_id"] = waypoints.apply(lambda x: x["milestone_id"] if x["percentage"] == 1 else None, axis=1) # 不在里程碑上的waypoint的milestone_id为None
    waypoints = waypoints[["waypoint_id", "milestone_id"]]

    wp = {
        "milestone_percentages": waypoint_milestone_percentages,
        "progressions": waypoint_progressions,
        "geodesic_distances": waypoint_geodesic_distances,
        "waypoint_network": waypoint_network,
        "waypoints": waypoints,
    }

    return wp


def is_wrapper_with_waypoints(trajectory):
    return trajectory.get("pydynwrap:with_waypoints", False)
