import pandas as pd
import networkx as nx

import pydynverse as pdv
from .simplify_networkx_network import simplify_networkx_network


def simplify_trajectory(trajectory, allow_self_loops=False):
    # TODO: 简化轨迹
    # 构造igraph图
    is_directed = trajectory["milestone_network"]["directed"].all()
    gr = nx.from_pandas_edgelist(
        trajectory["milestone_network"],
        source="from",
        target="to",
        edge_attr=True,
        create_using=nx.DiGraph if is_directed else nx.Graph
    )

    # 简化细胞
    edge_points = trajectory["progressions"]
    edge_points.rename(columns={"cell_id": "id"}, inplace=True)
    edge_points["id"] = edge_points["id"].apply(lambda x: f"SIMPLIFYCELL_{x}")

    # 核心操作：简化igraph网络结构
    out = simplify_networkx_network(
        gr,
        allow_duplicated_edges=False,
        allow_self_loops=allow_self_loops,
        force_keep=trajectory["divergence_regions"]["milestone_id"],
        edge_points=edge_points
    )  # 暂时只输出简化后的networkx图结构

    # 基于简化后的igraph图结构milestone相关数据结构
    gr = out["gr"]
    milestone_ids = gr.nodes
    milestone_network = pd.DataFrame(gr.edges(data=True), columns=["from", "to", "attributes"])
    milestone_network = pd.concat([milestone_network.drop(columns=['attributes']), milestone_network["attributes"].apply(pd.Series)], axis=1)
    milestone_network = milestone_network[["from", "to", "weight", "directed"]].rename(columns={"weight": "length"})

    edge_points = out["edge_points"]
    progressions = out["edge_points"][["id", "from", "to", "percentage"]].rename(columns={"id": "cell_id"})
    progressions["cell_id"] = progressions["cell_id"].apply(lambda x: x.replace("SIMPLIFYCELL_", ""))

    # new_trajectory = trajectory.copy()
    pdv.wrap.add_trajectory(
        dataset=trajectory,
        milestone_ids=milestone_ids,
        milestone_network=milestone_network,
        divergence_regions=trajectory["divergence_regions"],
        progressions=progressions,
        allow_self_loops=allow_self_loops,
    )

    # TODO: 带降维的细胞简化
    return trajectory
