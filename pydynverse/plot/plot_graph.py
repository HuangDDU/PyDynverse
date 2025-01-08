import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import scanpy as sc

from ..wrap import calculate_trajectory_dimred
from .add_cell_coloring import add_cell_coloring
from .add_milestone_coloring import add_milestone_coloring


def plot_graph(trajectory, color_cells="milestone", color_milestones=None):
    dimred_traj = calculate_trajectory_dimred(trajectory)

    # 添加里程碑milestone颜色
    milestones = add_milestone_coloring(trajectory["milestone_ids"], color_milestones)

    #  添加细胞颜色
    cell_positions = dimred_traj["cell_positions"]
    cell_coloring_output = add_cell_coloring(
        cell_positions=cell_positions,
        trajectory=trajectory,
        color_cells=color_cells,
        milestones=milestones
    )

    milestone_positions = dimred_traj["milestone_positions"]

    # 绘制细胞
    adata = trajectory["adata"]
    cell_id_list = list(adata.obs.index)
    adata.obsm["dimred"] = cell_positions.loc[cell_id_list].values
    color_cells = cell_coloring_output["color_cells"]
    if color_cells == "milestone":
        # 先绘制网络
        ax = plt.subplots(1, 1)[1]
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
        )
        # 再绘制细胞
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
