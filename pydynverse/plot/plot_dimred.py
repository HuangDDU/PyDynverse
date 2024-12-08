import pandas as pd

import scanpy as sc

from .dummy_proofing import check_milestones
from .add_milestone_coloring import add_milestone_coloring
from .add_cell_coloring import add_cell_coloring
from .project_waypoints import project_waypoints_coloured

from ..wrap.wrap_add_dimred import get_dimred

def plot_dimred(
        trajectory,
        color_cells="auto",
        dimred="dimred",
        plot_trajectory=True,
        plot_milestone_network=True,
        size_milestones=6,
        size_transitions=2,
        # 轨迹信息
        grouping=None,
        feature_oi=None,
        color_milestones=None,
        groups=None,
        milestones=None,
        milestone_percentages=None,
        pseudotime=None,
        expression_source="expression",
):
    # TODO: 提取里程碑
    milestones = check_milestones(trajectory, milestones=milestones, milestone_percentages=milestone_percentages)
    dimred, dimred_extra = get_dimred(
        dataset=trajectory,
        dimred=dimred,
        expression_source=expression_source,
        return_other_dimreds=True
    )  # 提取降维结果, 这里与R不同, 同时返回dimred, dimred_extra

    # TODO: 提取细胞位置, 基于降维结果与里程碑百分比
    cell_positions = dimred

    # TODO:添加里程碑颜色
    if (plot_milestone_network or plot_trajectory) and color_cells == "milestone":
        if not "color" in milestones:
            milestones = add_milestone_coloring(milestones, color_milestones=color_milestones)

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
    adata = trajectory["adata"]
    adata.obs["grouping"] = grouping
    adata.uns["grouping_colors"] = cell_coloring_output["color_dict"]
    adata.obsm["dimred"] = dimred.values
    ax = sc.pl.embedding(adata, basis="dimred", color="grouping", show=False) # 后需要在这个画板上继续添加内容

    # TODO: 绘制里程碑，箭头等

    # 绘制轨迹
    if plot_trajectory:
        edge_positions = None
        dimred_segment_progressions = dimred_extra["dimred_segment_progressions"]
        dimred_segment_points = dimred_extra["dimred_segment_points"]
        if not (dimred_segment_progressions is None) and not (dimred_segment_points is None):
            edge_positions = pd.concat([dimred_segment_progressions, dimred_segment_points], axis=1)
            # TODO: 需要了解WayPoint概念,真正的投影
            milestone_positions = edge_positions[edge_positions["percentage"]==0]
            ax.scatter(edge_positions[0], edge_positions[1], c="black", s=size_transitions)
            ax.scatter(milestone_positions[0], milestone_positions[1], c="black", s=size_milestones)
            
            # waypoint_projection = project_waypoints_coloured()
            # wp_segments = waypoint_projection # 到轨迹上的投影
            # milestone_positions = wp_segments # 过滤

            # 此处需要里程碑结构，

            # color = None
            # if grouping:
            #     color = grouping
            # elif feature_oi:
            #     color = feature_oi
            # 提取AnnData方便绘制细胞降维图
            # sc.pl.embedding(adata, basis=dimred)
