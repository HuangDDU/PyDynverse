import pandas as pd

import scanpy as sc

from .dummy_proofing import check_milestones
from .add_milestone_coloring import add_milestone_coloring
from .add_cell_coloring import add_cell_coloring
from .project_waypoints import project_waypoints_coloured

from ..wrap.wrap_add_dimred import get_dimred, is_wrapper_with_dimred
from ..dimred import dimred_mds, dimred_tsne, dimred_iosmap


def plot_dimred(
        trajectory,
        color_cells="auto",
        dimred=None,
        plot_trajectory=True,
        plot_milestone_network=True,
        size_milestones=30,
        size_transitions=2,
        # 轨迹信息
        grouping=None,
        feature_oi=None,
        color_milestones=None,
        groups=None,
        milestones=None,
        milestone_percentages=None,
        pseudotime=None,
        expression_source="counts",
        arrow=None,
        # 投影轨迹的相关参数
        waypoints=None,
        trajectory_projection_sd=None,
        color_trajectory=None
):
    # 提取里程碑
    milestones = check_milestones(trajectory, milestones=milestones, milestone_percentages=milestone_percentages)

    # 提取降维结果
    if (dimred is None) and (not is_wrapper_with_dimred(trajectory)):
        # Note: test dimred method for simulation data
        dimred = dimred_mds
        # dimred = dimred_tsne
        # dimred = dimred_iosmap
    dimred, dimred_extra = get_dimred(
        dataset=trajectory,
        dimred=dimred,
        expression_source=expression_source,
        return_other_dimreds=True
    )  # 提取降维结果, 这里与R不同, 同时返回dimred, dimred_extra

    #  提取细胞位置, 基于降维结果与里程碑百分比
    cell_positions = dimred
    cell_positions["cell_id"] = cell_positions.index
    # merge链接操作降维结果与里程碑百分比调整细胞位置
    tmp_milestone_percentages = trajectory["milestone_percentages"]
    idx = tmp_milestone_percentages.groupby("cell_id")["percentage"].idxmax()
    tmp_milestone_percentages = tmp_milestone_percentages.loc[idx, ["cell_id", "milestone_id"]]
    cell_positions = pd.merge(cell_positions, tmp_milestone_percentages, on="cell_id")

    # 添加里程碑颜色
    if (plot_milestone_network or plot_trajectory) and color_cells == "milestone":
        milestones = add_milestone_coloring(trajectory["milestone_ids"], color_milestones=color_milestones)

    # 最主要部分, 添加细胞颜色
    cell_coloring_output = add_cell_coloring(
        cell_positions=cell_positions,
        trajectory=trajectory,
        color_cells=color_cells,
        grouping=grouping,
        groups=groups,
        feature_oi=feature_oi,
        expression_source=expression_source,
        pseudotime=pseudotime,
        color_milestones=color_milestones,
        milestones=milestones,
        milestone_percentages=milestone_percentages
    )  # NOTE:只是为了获得一个颜色字典，这么费劲, 再CFE可以优化

    # TODO: 没有颜色则计算细胞密度

    # 绘制细胞
    color_cells = cell_coloring_output["color_cells"]
    adata = trajectory["adata"]
    adata.obsm["dimred"] = dimred.values
    if color_cells == "grouping":
        adata.obs["grouping"] = grouping
        adata.uns["grouping_colors"] = cell_coloring_output["color_dict"]
        ax = sc.pl.embedding(adata, basis="dimred", color="grouping", show=False)  # 后需要在这个画板上继续添加内容
    elif color_cells == "milestone":
        cell_id_list = adata.obs.index
        color_list = cell_coloring_output["color_scale"][cell_id_list].to_list()
        adata.obs["cell_id"] = pd.Categorical(cell_id_list, categories=cell_id_list)
        adata.uns["cell_id_colors"] = color_list
        ax = sc.pl.embedding(adata, basis="dimred", color="cell_id", show=False, legend_loc=None)

    # TODO: 绘制里程碑，箭头等
    # 绘制轨迹
    if plot_trajectory:
        dimred_segment_progressions = dimred_extra.get("dimred_segment_progressions", None)
        dimred_segment_points = dimred_extra.get("dimred_segment_points", None)
        if not (dimred_segment_progressions is None) and not (dimred_segment_points is None):
            edge_positions = pd.concat([dimred_segment_progressions, dimred_segment_points], axis=1)
            # # TODO: 需要了解WayPoint概念,真正的投影
            # milestone_positions = edge_positions[edge_positions["percentage"] == 0]
            # ax.scatter(edge_positions[0], edge_positions[1], c="black", s=size_transitions)
            # ax.scatter(milestone_positions[0], milestone_positions[1], c="black", s=size_milestones)
        else:
            edge_positions = None

        waypoint_projection = project_waypoints_coloured(
            trajectory=trajectory,
            cell_positions=cell_positions,
            waypoints=waypoints,
            trajectory_projection_sd=trajectory_projection_sd,
            color_trajectory=color_trajectory,
            edge_positions=edge_positions
        )
        wp_segments = waypoint_projection["segments"]  # 到轨迹上的投影
        milestone_positions = wp_segments[wp_segments["milestone_id"].apply(lambda x: not x is None)]  # 过滤只保留里程碑节点上的waypoint点, 就是milestone节点

        # 绘制waypoint构成的轨迹曲线
        ax.scatter(milestone_positions["comp_1"], milestone_positions["comp_2"], c="black", s=size_milestones)  # 绘制里程碑附近的点
        # 分段绘制waypoint连成的曲线
        for g in wp_segments["group"].unique():
            wp_segments_g = wp_segments[wp_segments["group"] == g]
            ax.plot(wp_segments_g["comp_1"], wp_segments_g["comp_2"], c="black", linewidth=size_transitions)

        # 绘制轨迹上的箭头
        # TODO: 箭头大小调整
        if trajectory["milestone_network"]["directed"].any():
            def get_arrow_df(group):
                group = group.sort_values(by="percentage")
                start = group.iloc[0]
                end = group.iloc[-1]
                s = pd.Series({
                    "x": start["comp_1"],
                    "y": start["comp_2"],
                    "dx": end["comp_1"]-start["comp_1"],
                    "dy": end["comp_2"]-start["comp_2"]}
                )
                return s
            arrow_df = wp_segments[wp_segments["arrow"]].groupby("group").apply(get_arrow_df)
            ax.quiver(arrow_df["x"], arrow_df["y"], arrow_df["dx"], arrow_df["dy"])
            if color_trajectory is None:
                # TODO: 轨迹添加颜色
                pass
            else:
                pass

        # TODO: 绘制waypoint分段

        # TODO: 添加里程碑标签
