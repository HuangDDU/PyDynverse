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
    pos = nx.spring_layout(gr)
    layout = pd.DataFrame(pos).T  # TODO: 布局可能还需要Z归一化

    # 样本的降维投影策略
    def mix_dimred(milid, milpct):
        return layout.loc[milid].apply(lambda x: (x.array * milpct.array).sum())  # 降维的每个维度按比例确定值

    # 对样本细胞降维
    cell_positions = milestone_percentages.groupby("cell_id").apply(lambda x: mix_dimred(x["milestone_id"], x["percentage"]))

    milestone_positions = layout

    edge_positions = None  # TODO: 这里暂时不需要手动绘制边，NetworkX会自动绘制边

    # TODO: 提取发散区域的线或多边形区域
    if (divergence_regions is not None) and (divergence_regions.shape[0] > 0):
        pass
    else:
        pass

    divergence_edge_positions = None
    divergence_polygon_positions = None

    return {
        "milestone_positions": milestone_positions,
        "edge_positions": edge_positions,  # None
        "cell_positions": cell_positions,
        "divergence_edge_positions": divergence_edge_positions,  # None
        "divergence_polygon_positions": divergence_polygon_positions,  # None
        # 额外添加
        "gr": gr,
        "pos": pos,
    }


def get_divergence_triangles(divergence_regions):
    return None
