import numpy as np
import pandas as pd
import networkx as nx
# from .wrap_add_trajectory import is_wrapper_with_trajectory
# from .wrap_add_waypoints import is_wrapper_with_waypoints
from ..util import calculate_distance


def calculate_geodesic_distances(
    trajectory,
    waypoint_cells=None,
    waypoint_milestone_percentages=None,
    directed=False
):
    # assert is_wrapper_with_trajectory(trajectory)

    # if (waypoint_cells is None) and is_wrapper_with_waypoints(trajectory):
    #     waypoint_cells <- trajectory["waypoint_cells"]
    return calculate_geodesic_distances_(
        cell_ids=trajectory["cell_ids"],
        milestone_ids=trajectory["milestone_ids"],
        milestone_network=trajectory["milestone_network"],
        milestone_percentages=trajectory["milestone_percentages"],
        divergence_regions=trajectory["divergence_regions"],
        waypoint_cells=waypoint_cells,
        waypoint_milestone_percentages=waypoint_milestone_percentages,
        directed=directed
    )


def calculate_geodesic_distances_(
    cell_ids,
    milestone_ids,
    milestone_network,
    milestone_percentages,
    divergence_regions,
    waypoint_cells=None,
    waypoint_milestone_percentages=None,
    directed=None
):
    # 整体思路： 分别在每个发散区域内部计算目标点的全路径， 在整体图上合并计算距离
    # NOTE: 0. 预处理
    # 传参进来的cell_ids可能是部分,这里提取轨迹的全部cell_ids
    cell_ids_trajectory = list(set(milestone_percentages["cell_id"]))

    # 提取waypoint和milestone percentages
    if not waypoint_cells is None:
        waypoint_ids = waypoint_cells
    elif waypoint_milestone_percentages is None:
        waypoint_ids = cell_ids_trajectory
    else:
        waypoint_ids = []

    if not waypoint_milestone_percentages is None:
        waypoint_ids = waypoint_ids + waypoint_milestone_percentages["waypoint_id"].unique().tolist()
        milestone_percentages = pd.concat([
            milestone_percentages,
            waypoint_milestone_percentages.rename(columns={"waypoint_id": "cell_id"})
        ])  # 叠加milestone和waypoint

    if divergence_regions is None:
        divergence_regions = pd.DataFrame(columns=["divergence_id", "milestone_id", "is_start"])

    # 重命名milestones_id来避免冲突, 添加前缀并在milestone_id相关内容里修改
    def milestone_trafo_fun(x): return f"MILESTONE_{x}"
    milestone_network["from"] = milestone_network["from"].apply(milestone_trafo_fun)
    milestone_network["to"] = milestone_network["to"].apply(milestone_trafo_fun)
    milestone_ids = list(map(milestone_trafo_fun, milestone_ids))
    milestone_percentages["milestone_id"] = milestone_percentages["milestone_id"].apply(milestone_trafo_fun)
    divergence_regions["milestone_id"] = divergence_regions["milestone_id"].apply(milestone_trafo_fun)

    # 添加extra发散区域, 正常的边也被当作发散区域
    extra_divergences = milestone_network.copy()
    extra_divergences = extra_divergences[~(extra_divergences["from"] == extra_divergences["to"])]  # 过滤自环
    # in_divergence判断当前边是否在已有的发散区域内，标准的延迟承诺区域
    divergence_regions_set_list = divergence_regions.groupby("divergence_id")["milestone_id"].apply(set).tolist()

    def is_milestone_in_divergence(milestone_set, divergence_regions_set_list):
        for divergence_regions_set in divergence_regions_set_list:
            if milestone_set.issubset(divergence_regions_set):
                return True
        return False
    extra_divergences["in_divergence"] = extra_divergences.apply(lambda x: is_milestone_in_divergence({x["from"], x["to"]}, divergence_regions_set_list), axis=1)
    extra_divergences = extra_divergences[~extra_divergences["in_divergence"]]  # 只保留新的发散区域
    extra_divergences["divergence_id"] = extra_divergences.apply(lambda x: f"{x['from']}__{x['to']}", axis=1)
    extra_divergences = pd.concat([
        # 添加新的milestone_id和is_start列
        extra_divergences.assign(milestone_id=extra_divergences["from"], is_start=True),
        extra_divergences.assign(milestone_id=extra_divergences["to"], is_start=False)
    ])[["divergence_id", "milestone_id", "is_start"]]

    # 合并发散区域
    divergence_regions = pd.concat([divergence_regions, extra_divergences]).reset_index(drop=True)
    divergence_regions["is_start"] = divergence_regions["is_start"].astype(bool) # 确保is_start为bool类型
    divergence_ids = divergence_regions["divergence_id"].unique()

    # 准备使用NetworkX, 构造相关数据, 从DataFrame开始构造
    milestone_graph = nx.from_pandas_edgelist(milestone_network, source="from", target="to", edge_attr="length")

    # NOTE: 1. 分别计算
    # 计算发散区域内部细胞间距离
    def calc_divergence_inner_distance_df(did):
        dir = divergence_regions[divergence_regions["divergence_id"] == did]
        mid = dir[dir["is_start"]]["milestone_id"].tolist()  # 该区域起点milestone_id
        tent = dir["milestone_id"].tolist()  # 该区域所有milestone_id
        tent_distances = pd.DataFrame(index=mid, columns=tent, data=np.zeros((len(mid), len(tent))))  # 区域内起点到所有milestone的距离
        # 从图中提取对应的边
        for i in mid:
            for j in tent:
                if i == j:
                    tent_distances.loc[i, j] = 0
                else:
                    tent_distances.loc[i, j] = milestone_graph.edges[(i, j)]["length"]
        # 此处复用来找相关的点的cell_id
        relevant_pct_cell_id_list = milestone_percentages.groupby("cell_id")["milestone_id"].apply(lambda x: is_milestone_in_divergence(set(x), [set(tent)]))
        relevant_pct_cell_id_list = relevant_pct_cell_id_list[relevant_pct_cell_id_list].index.to_list()
        relevant_pct = milestone_percentages[milestone_percentages["cell_id"].apply(lambda x: x in relevant_pct_cell_id_list)]
        if relevant_pct.shape[0] <= 1:
            return None

        scaled_dists = relevant_pct.copy()
        scaled_dists["dist"] = scaled_dists.apply(lambda x: x["percentage"]*tent_distances.loc[mid, x["milestone_id"]], axis=1)
        tent_distances_long = tent_distances.melt(var_name="from", value_name="length")  # 宽数据转化为长数据
        tent_distances_long["to"] = tent_distances_long["from"]

        pct_mat = pd.concat([
            scaled_dists[["cell_id", "milestone_id", "dist"]].rename(columns={"cell_id": "from", "milestone_id": "to", "dist": "length"}),
            tent_distances_long
        ]).pivot(index="from", columns="to", values="length").fillna(0)  # (n_cell+n_milestone+n_waypoint)*n_milestone, 长数据转宽数据, 索引名为from

        wp_cells = list(set(pct_mat.index) & set(waypoint_ids))

        if directed:
            # TODO: 暂时不管有向图
            pass

        distances = calculate_distance(pct_mat, pct_mat.loc[wp_cells+tent], method="manhattan")
        # TODO: 这里calculate_distance可以直接用sklearn的pairwise_distances代替
        distances = distances.reset_index().melt(id_vars="from", var_name="to",  value_name="length")  # 宽数据转化为长数据
        distances = distances[~(distances["from"] == distances["to"])]
        return distances

    cell_in_tent_distances = pd.concat([calc_divergence_inner_distance_df(did) for did in divergence_ids])

    if directed:
        # TODO: 暂时不管有向图
        pass

    # NOTE: 2. 合并计算
    # 合并两个图到一张图上
    graph = pd.concat([milestone_network, cell_in_tent_distances]).groupby(["from", "to"]).agg({"length": "min"}).reset_index()  # 合并后提取最短边，目前没什么效果，可能对环图有用
    graph = nx.from_pandas_edgelist(graph, source="from", target="to", edge_attr="length")
    # 选择后续有向图的最短距离模式
    if directed or directed == "forward":
        mode = "out"
    elif directed == "reverse":
        mode = "in"
    else:
        mode = "all"

    # 调用dijkstra计算, 计算waypoint和cell之间的
    out = pd.DataFrame(np.zeros((len(waypoint_ids), len(cell_ids))), index=waypoint_ids, columns=cell_ids)
    for source in waypoint_ids:
        length_dict = nx.single_source_dijkstra_path_length(graph, source=source, weight="length")
        for target in cell_ids:
            if target in length_dict:
                out.loc[source, target] = length_dict[target]
            else:
                out.loc[source, target] = np.inf

    # TODO: 过滤一些细胞
    cell_ids_filtered = []
    if len(cell_ids_filtered) > 0:
        pass

    return out.loc[waypoint_ids, cell_ids]


def calculate_geodesic_distances_me():
    # TDOO: 整体构造NetworkX图，计算最短距离
    pass
