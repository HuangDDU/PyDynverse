import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import networkx as nx
import scanpy as sc

from ..wrap import calculate_trajectory_dimred
from .add_cell_coloring import add_cell_coloring
from .add_milestone_coloring import add_milestone_coloring


def plot_graph(trajectory, color_cells="auto", color_milestones=None, grouping=None):
    dimred_traj = calculate_trajectory_dimred(trajectory)

    # 添加里程碑milestone颜色
    milestones = add_milestone_coloring(trajectory["milestone_ids"], color_milestones)

    #  添加细胞颜色
    cell_positions = dimred_traj["cell_positions"]
    cell_coloring_output = add_cell_coloring(
        cell_positions=cell_positions,
        trajectory=trajectory,
        color_cells=color_cells,
        grouping=grouping,
        milestones=milestones
    )

    milestone_positions = dimred_traj["milestone_positions"]

    ax = plt.subplots(1, 1)[1]
    # 绘制分支承诺区域
    # 发散区域灰色背景
    dpp = dimred_traj["divergence_polygon_positions"]
    for triangle_id in dpp["triangle_id"].unique():
        polygon_vertices = dpp[dpp["triangle_id"] == triangle_id][["comp_1", "comp_2"]].values  # 提取边界点
        polygon = patches.Polygon(polygon_vertices, closed=True, fill=True, color="lightblue", alpha=0.5)
        ax.add_patch(polygon)
    # 发散区域灰色虚线
    dep = dimred_traj["divergence_edge_positions"]
    x_edges = dep[["comp_1_from", "comp_1_to"]].T.values  # 2*n
    y_edges = dep[["comp_2_from", "comp_2_to"]].T.values  # 2*n
    ax.plot(x_edges, y_edges, color="gray", linestyle="--", linewidth=5)

    # 绘制网络
    G = dimred_traj["gr"]
    pos = dimred_traj["pos"]
    milestone_size = 300
    nx.draw(
        G,
        pos,
        with_labels=True,
        node_color=[milestones.loc[node, "color"] for node in G.nodes],
        node_size=milestone_size,
        width=10,
        edge_color="gray",
    )  # TODO: 边的样式更改。

    # 再绘制细胞
    adata = trajectory["adata"]
    cell_id_list = list(adata.obs.index)
    adata.obsm["dimred"] = cell_positions.loc[cell_id_list].values
    color_cells = cell_coloring_output["color_cells"]
    if color_cells == "grouping":
        adata.obs["grouping"] = grouping
        adata.uns["grouping_colors"] = cell_coloring_output["color_dict"]
        ax = sc.pl.embedding(
            adata,
            basis="dimred",
            color="grouping",
            show=False,
            title="",
            ax=ax
        )  # 后需要在这个画板上继续添加内容
    elif color_cells == "milestone":
        color_list = cell_coloring_output["color_scale"][cell_id_list].to_list()
        adata.obs["cell_id"] = pd.Categorical(cell_id_list, categories=cell_id_list)
        adata.uns["cell_id_colors"] = color_list
        cell_size = milestone_size
        ax = sc.pl.embedding(
            adata,
            basis="dimred",
            color="cell_id",
            show=False,
            legend_loc=None,
            size=cell_size,
            title="",
            ax=ax
        )
