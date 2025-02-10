import itertools
import pandas as pd
import networkx as nx


def calculate_trajectory_dimred(trajectory, adjust_weights=False):
    # 一堆assert判断参数有效性，先跳过

    # 从trajectory提取数据
    cell_ids = trajectory["cell_ids"]
    milestone_ids = trajectory["milestone_ids"]
    num_milestones = len(milestone_ids)
    milestone_network = trajectory["milestone_network"]
    milestone_percentages = trajectory["milestone_percentages"]
    is_directed = trajectory["milestone_network"]["directed"].any()

    structure = milestone_network.copy()

    # 在发散区域的每两个终端节点对之间添加不可见的边
    divergence_regions = trajectory.get("divergence_regions", None)
    if (divergence_regions is not None) and (divergence_regions.shape[0] > 0):
        # TODO:
        pass
        divergence_edges = get_divergence_triangles(divergence_regions)
        # structure = pd.concat([structure, divergence_edges])
        structure = pd.concat([structure])

    if adjust_weights:
        # TODO:
        # 权重调整
        pass

    # 添加边权重
    structure["weight"] = structure["length"]  # TODO: 暂时直接把边长作为边权，后续可能需要调整

    # 使用networkx包进行milestone network的节点布局
    gr = nx.from_pandas_edgelist(structure, source="from", target="to", edge_attr=True, create_using=nx.DiGraph if is_directed else nx.Graph)
    # pos = nx.kamada_kawai_layout(gr)
    # pos = nx.spring_layout(gr, seed=0)  # dict node:position
    pos = nx.nx_agraph.graphviz_layout(gr, prog="dot")  # graphviz的dot有向无环图布局更加合适
    layout = pd.DataFrame(pos).T  # dataframe
    # 是否Z归一化差距不大
    # layout = pd.DataFrame(StandardScaler().fit_transform(layout), columns=layout.columns, index=layout.index)  # 布局Z归一化
    # pos = dict(zip(layout.index, layout.values))

    # 样本的降维投影策略
    def mix_dimred(milid, milpct):
        return layout.loc[milid].apply(lambda x: (x.array * milpct.array).sum())  # 降维的每个维度按比例确定值

    # 对样本细胞降维
    cell_positions = milestone_percentages.groupby("cell_id").apply(lambda x: mix_dimred(x["milestone_id"], x["percentage"]))

    milestone_positions = layout

    edge_positions = None  # NOTE: 这里暂时不需要手动绘制边，NetworkX会自动绘制边

    # 提取发散区域的线或多边形区域
    if (divergence_regions is not None) and (divergence_regions.shape[0] > 0):
        triags = get_divergence_triangles(divergence_regions)
        # 边坐标稍后绘制虚线
        divergence_edge_positions = triags.rename(columns={"node1": "from", "node2": "to"})
        divergence_edge_positions[["comp_1_from", "comp_2_from"]] = divergence_edge_positions["from"].apply(lambda x: milestone_positions.loc[x])
        divergence_edge_positions[["comp_1_to", "comp_2_to"]] = divergence_edge_positions["to"].apply(lambda x: milestone_positions.loc[x])
        # 多边形区域稍后绘制阴影
        divergence_polygon_positions = triags.copy()
        divergence_polygon_positions["triangle_id"] = [f"triangle_{i}" for i in range(triags.shape[0])]
        divergence_polygon_positions = divergence_polygon_positions.melt(
            id_vars=["triangle_id"],
            value_vars=["start", "node1", "node2"],
            var_name="triangle_part",
            value_name="milestone_id"
        )  # 沿着triangle_id展开宽数据为长数据
        divergence_polygon_positions[["comp_1", "comp_2"]] = divergence_polygon_positions["milestone_id"].apply(lambda x: milestone_positions.loc[x])
    else:
        divergence_edge_positions = pd.DataFrame(columns=["divergence_id", "start", "from", "to", "comp_1_from", "comp_2_from", "comp_1_to", "comp_2_to"])
        divergence_polygon_positions = pd.DataFrame(columns=["triangle_id", "triangle_part", "milestone_id", "comp_1", "comp_2"])

    return {
        "milestone_positions": milestone_positions,
        "edge_positions": edge_positions,  # None
        "cell_positions": cell_positions,
        "divergence_edge_positions": divergence_edge_positions,
        "divergence_polygon_positions": divergence_polygon_positions,
        # 额外添加
        "gr": gr,
        "pos": pos,
    }


def get_divergence_triangles(divergence_regions):
    # TODO: 在CellFateExplorer中可以替换成循环实现
    def get_divergence_edge_df(did):
        rel_did = divergence_regions[divergence_regions["divergence_id"] == did]

        fr = rel_did[rel_did["is_start"]]["milestone_id"].tolist()[0]  # 只有一个
        tos = rel_did[~rel_did["is_start"]]["milestone_id"].tolist()

        de_df = pd.DataFrame(itertools.product(tos, tos), columns=["node1", "node2"])
        de_df = de_df[de_df["node1"] < de_df["node2"]]
        de_df["divergence_id"] = did
        de_df["start"] = fr

        return de_df

    triangles = pd.concat([get_divergence_edge_df(did) for did in divergence_regions["divergence_id"].unique()])

    return triangles[["divergence_id", "start", "node1", "node2"]]
