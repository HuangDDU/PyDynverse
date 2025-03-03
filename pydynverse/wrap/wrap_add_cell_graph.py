import numpy as np
import pandas as pd
import networkx as nx
import igraph as ig

from .wrap_add_trajectory import add_trajectory
from .simplify_trajectory import simplify_trajectory


def add_cell_graph(
        dataset,
        cell_graph: pd.DataFrame,
        to_keep: pd.Series | dict = None,
        milestone_prefix: str = "milestone_",
        backend: str = "networkx"
):
    if not "length" in cell_graph.columns:
        cell_graph["length"] = 1
    if not "directed" in cell_graph.columns:
        cell_graph["directed"] = False

    cell_ids = dataset["cell_ids"]
    is_directed = cell_graph["directed"].any()

    # 关键节点
    if to_keep is None:
        to_keep = pd.Series(True, index=dataset.cell_ids)
    elif type(to_keep) == dict:
        to_keep = pd.Series(to_keep)
    v_keeps = to_keep[to_keep].index.to_list()

    if backend.lower() == "networkx":
        # 构造图, networkX存取dataframe更加方便
        G = nx.from_pandas_edgelist(cell_graph, source="from", target="to", edge_attr=["length", "directed"], create_using=nx.DiGraph if is_directed else nx.Graph)  # 构造network graph对象

        # 初步简化该图
        # STEP 1: for each cell, find closest milestone
        # 步骤1： 对于每个细胞，找到最近的里程碑
        distance_df = pd.DataFrame(dict(nx.shortest_path_length(G.to_undirected(),  weight="length"))).loc[cell_ids, v_keeps]  # 计算距离时当作无向边考虑，相当于igraph的mode=all
        closest_trajpoint = distance_df.idxmin(axis=1)  # 每个i细胞对应的最近的关键节点

        # STEP 2: simplify backbone
        # 步骤2：简化骨架，诱导子图: 子图中两两顶点在原图中的边一定在子图中存在
        G = G.subgraph(v_keeps)
        milestone_ids = G.nodes

        # STEP 3: Calculate progressions of cell_ids to determine which nodes were on each path
        # 步骤3：计算细胞的progressions来决定点在哪个路径上
        milestone_network_proto = nx.to_pandas_edgelist(G, source="from", target="to")
        milestone_network_proto["path"] = milestone_network_proto.apply(lambda x: nx.shortest_path(G, source=x["from"], target=x["to"]), axis=1)
        # 计算关键节点的percentage
        progressions_v_keeps = milestone_network_proto\
            .explode("path")\
            .groupby("path")\
            .agg(lambda x: x.iloc[0]).reset_index()\
            .rename(columns={"path": "node"})[["from", "to", "length", "node"]]  # 保留关键节点所在的第一条边
        progressions_v_keeps["percentage"] = progressions_v_keeps.apply(lambda x: nx.shortest_path_length(G, source=x["from"], target=x["node"],  weight="length")/x["length"], axis=1)

        closest_trajpoint_df = pd.DataFrame()
        closest_trajpoint_df["node"] = closest_trajpoint
        closest_trajpoint_df["cell_id"] = cell_ids
        progressions = pd.merge(progressions_v_keeps, closest_trajpoint_df, on="node")  # 拼接实现映射所有细胞到最近的关键节点上
        progressions = progressions[["cell_id", "from", "to", "percentage"]]

        milestone_network = milestone_network_proto[["from", "to", "length", "directed"]]

        # 添加里程碑名称前缀
        milestone_ids = [f"{milestone_prefix}{milestone_id}" for milestone_id in milestone_ids]
        milestone_network[["from", "to"]] = milestone_prefix + milestone_network[["from", "to"]]
        progressions[["from", "to"]] = milestone_prefix + progressions[["from", "to"]]

    else:
        # TODO: igraph更快的实现
        # # make network graph object,
        # # 构造图, 后续要用到找最短距离的to指定节点，igraph更快
        # gr = ig.Graph.TupleList(cell_graph.values, edge_attrs=["length", "directed"], directed=is_directed)

        # # 初步简化该图
        # # STEP 1: for each cell, find closest milestone
        # # 步骤1： 对于每个细胞，找到最近的里程碑
        # dists = gr.distances(target=v_keeps, weights="length")
        # closest_trajpoint = [v_keeps[i] for i in np.argmin(dists, axis=1)]
        # # closest_trajpoint = pd.Series(data=closest_trajpoint, index=gr.vs["name"]) # 要用Series保存，graph中节点顺序可能与dataset节点顺序不一致
        # closest_trajpoint = dict(zip(gr.vs["name"], closest_trajpoint))  # 要用dict保存，graph中节点顺序可能与dataset节点顺序不一致

        # # STEP 2: simplify backbone
        # # 步骤2：简化骨架，诱导子图: 子图中两两顶点在原图中的边一定在子图中存在
        # gr = gr.induced_subgraph(v_keeps)
        # milestone_ids = gr.vs["name"]

        # # STEP 3: Calculate progressions of cell_ids to determine which nodes were on each path
        # # 步骤3：计算细胞的progression来决定点在哪个路径上
        # milestone_network_proto = pd.DataFrame()
        # # for each node, find an edge which contains the node and calculate its progression along that edge
        # # 对于每个节点，找到所在边并计算沿边的progression
        # progressions = milestone_network_proto  # v_keep的progressions
        # # pd.merge() # 拼接补充其他的progression
        milestone_network = None
        progressions = None

    trajectory = add_trajectory(
        dataset=dataset,
        milestone_network=milestone_network,
        divergence_regions=None,
        progressions=progressions,
    )
    
    simplified_trajectory = simplify_trajectory(trajectory) # 最后简化轨迹
    return simplified_trajectory
