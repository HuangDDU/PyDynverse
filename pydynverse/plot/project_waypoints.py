import numpy as np
import pandas as pd
from scipy.stats import norm

from ..wrap import select_waypoints
from .._logging import logger


def project_waypoints_coloured(
    trajectory,
    cell_positions,
    edge_positions=None,
    waypoints=None,
    trajectory_projection_sd=None,
    color_trajectory=None
):
    if waypoints is None:
        # TODO: 真实数据集上，这里示例数据集在Slingshot运行的太慢大约40s，而在PAGA运行1s，需要分析原因并加速
        logger.debug("selecting waypoints")
        waypoints = select_waypoints(trajectory)
        logger.debug(f"selecting waypoints shape is {waypoints['geodesic_distances'].shape}, finished!")

    if trajectory_projection_sd is None:
        trajectory_projection_sd = sum(trajectory["milestone_network"]["length"])*0.05

    wps = waypoints
    wps["waypoint_network"] = wps["waypoint_network"].rename({"from_milestone_id": "milestone_id_from", "to_milestone_id": "milestone_id_to"})  # 简单换列名, 之后再也没出现过,没啥用

    # 计算waypoint低维坐标
    if not edge_positions is None:
        # TODO: 给定边坐标, 暂时还是用与else逻辑一致， 后续修改
        # 使用核和几何距离计算wps的低维降维
        # 计算权重
        weights = norm.pdf(wps["geodesic_distances"], scale=trajectory_projection_sd)  # 高斯核
        weights /= weights.sum(axis=1, keepdims=True)  # 行归一化权重
        # 获取细胞位置
        positions = cell_positions[["cell_id", "comp_1", "comp_2"]].set_index("cell_id")
        positions = positions.loc[wps["geodesic_distances"].columns]
        # 计算降维结果
        result = np.dot(weights, positions)
        result_df = pd.DataFrame(result, columns=["comp_1", "comp_2"])
        result_df["waypoint_id"] = wps["geodesic_distances"].index
        # 合并结果
        waypoint_positions = pd.merge(result_df, wps["waypoints"], on="waypoint_id")
    else:
        # 使用核和几何距离计算wps的低维降维
        # 计算权重
        weights = norm.pdf(wps["geodesic_distances"], scale=trajectory_projection_sd)  # 高斯核
        weights /= weights.sum(axis=1, keepdims=True)  # 行归一化权重
        # 获取细胞位置
        positions = cell_positions[["cell_id", "comp_1", "comp_2"]].set_index("cell_id")
        positions = positions.loc[wps["geodesic_distances"].columns]
        # 计算降维结果
        result = np.dot(weights, positions)
        result_df = pd.DataFrame(result, columns=["comp_1", "comp_2"])
        result_df["waypoint_id"] = wps["geodesic_distances"].index
        # 合并结果
        waypoint_positions = pd.merge(result_df, wps["waypoints"], on="waypoint_id")

    if color_trajectory == "nearest":
        # TODO： 依据最近的细胞给定颜色
        pass

    # 构造最后输出
    segments = pd.merge(waypoint_positions, wps["progressions"], on="waypoint_id")
    segments["group"] = segments.apply(lambda x: f"{x['from']}---{x['to']}", axis=1)

    def calculate_closest_and_arrow(group):
        # 选择一段milestone_network边上的重点绘制箭头
        closest_index = (group["percentage"] - 0.5).abs().idxmin()  # 找到percentage最接近0.5的索引
        group["arrow"] = (group.index == closest_index) | (group.index == closest_index + 1)  # 创建arrow列
        return group
    segments = segments.groupby("group").apply(calculate_closest_and_arrow).reset_index(drop=True)

    waypoint_projection = {"segments": segments}

    return waypoint_projection
