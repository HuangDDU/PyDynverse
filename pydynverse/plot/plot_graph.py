import matplotlib.pyplot as plt
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
    color_cells = cell_coloring_output["color_cells"]
    adata = trajectory["adata"]
    adata.obsm["dimred"] = cell_positions.values
    if color_cells == "milestone":
        # 先绘制网络
        ax = plt.subplots(1,1)[1]
        G = dimred_traj["gr"]
        pos = dimred_traj["pos"]
        nx.draw(G, pos, with_labels=True, node_color = [milestones.loc[node, "color"] for node in G.nodes], edge_color="gray", width=10)
        # 再绘制细胞
        color_list = cell_coloring_output["color_scale"][adata.obs.index].to_list()
        adata.obs["cell_id"] = adata.obs.index.astype("category")
        adata.uns["cell_id_colors"] = dict(zip(adata.obs.index, color_list))
        ax = sc.pl.embedding(adata, basis="dimred", color="cell_id", show=False, legend_loc=None, size=1000, ax=ax)
        
