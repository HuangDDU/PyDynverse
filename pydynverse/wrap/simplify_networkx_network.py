import pandas as pd
import networkx as nx


def simplify_networkx_network(
    gr: nx.Graph,
    allow_duplicated_edges=True,
    allow_self_loops=True,
    force_keep=[],
    edge_points=None
):
    # 基于NetworkX简化milestone network

    # TODO: 无向图的边可以互相反转
    if edge_points is not None:
        edge_points = edge_points

    # 重新命名确保不混乱
    gr = nx.relabel_nodes(gr, {name: f"#M#{name}" for name in gr.nodes})
    force_keep = [f"#M#{name}" for name in force_keep]

    # 边权没有则手动添加为1
    attribute_keys = list(gr.edges(data=True))[0][2].keys()  # 提取第一个边的属性键名， 所有边的属性键名一致
    if "weight" not in attribute_keys:
        for u, v in gr.edges():
            gr[u][v]["weight"] = 1

    # 边是否有方向未指定，则手动添加与gr保持一致
    is_directed = gr.is_directed()
    if "directed" not in attribute_keys:
        for u, v in gr.edges():
            gr[u][v]["directed"] = is_directed

    # 分别处理每个连通分量
    if is_directed:
        connected_component_list = list(nx.weakly_connected_components(gr))  # 有向图提取弱连通分量，连接即可
    else:
        connected_component_list = list(nx.connected_components(gr))
    simplified_graphs = []
    for connected_component in connected_component_list:
        subgr = gr.subgraph(connected_component).copy()
        simplified_subgraph = simplify_subgraph(subgr, is_directed, force_keep)
        simplified_graphs.append(simplified_subgraph)

    # 合并图，这里直接合并
    out_gr = nx.compose_all([i["subgr"] for i in simplified_graphs])
    out_gr = nx.relabel_nodes(out_gr, {name: name[3:] for name in out_gr.nodes})  # 名字改回去

    if edge_points is None:

        return out_gr  # 暂时只输出简化后的networkx图结构
    else:
        sep = None
        return {
            "gr": out_gr,
            "sep": sep
        }


def simplify_subgraph(subgr, is_directed, force_keep):
    # NOTE: 这里是简化的核心函数
    node_list = list(subgr.nodes)
    keep_v = simplify_determine_nodes_to_keep(subgr, is_directed, force_keep)  # 决定保留哪些节点True， 过滤哪些节点False
    sub_edge_points = None  # TODO:
    keep_v = keep_v  # TODO:
    num_vs = len(subgr.nodes)
    neighs = simplify_get_neighbours(subgr, is_directed)
    to_process = [not i for i in keep_v]
    for v_rem in range(num_vs):
        if to_process[v_rem]:
            to_process[v_rem] = False
            # 入度前驱节点、边处理
            i = simplify_get_i(neighs, v_rem, is_directed)  # 前驱节点
            i_prev = v_rem
            left_path = [{"from": i, "to": i_prev, "weight": simplify_get_edge(subgr, i, i_prev)["weight"]}]
            while to_process[i]:
                # 前驱节点仍然需要删除
                to_process[i] = False
                tmp = i
                i = simplify_get_next(neighs, i, is_directed, left=True, prev=i_prev)
                i_prev = tmp
                # left_path.append({"from": i, "to": i_prev, "weight": simplify_get_edge(subgr, i, i_prev)["weight"]})
                left_path.append({"from": i, "to": i_prev, "weight": subgr.edges[(node_list[i], node_list[i_prev])]["weight"]})

            # 出度的后继节点边处理
            j = simplify_get_j(neighs, v_rem, is_directed)  # 后继节点
            j_prev = v_rem
            right_path = [{"from": j_prev, "to": j, "weight": simplify_get_edge(subgr, j_prev, j)["weight"]}]
            while to_process[j]:
                # 后继节点仍然需要删除
                to_process[j] = False
                tmp = j
                j = simplify_get_next(neighs, j, is_directed, left=False, prev=j_prev)
                j_prev = tmp
                # right_path.append({"from": j_prev, "to": j, "weight": simplify_get_edge(subgr, j_prev, j)["weight"]})
                right_path.append({"from": j_prev, "to": j, "weight": subgr.edges[(node_list[j_prev], node_list[j])]["weight"]})

            # 拼接后，节点序号转换为节点名字
            left_path = pd.DataFrame(left_path)
            left_path["from"] = [node_list[i] for i in left_path["from"]]
            left_path["to"] = [node_list[i] for i in left_path["to"]]
            right_path = pd.DataFrame(right_path)
            right_path["from"] = [node_list[i] for i in right_path["from"]]
            right_path["to"] = [node_list[i] for i in right_path["to"]]

            if i == j:
                # TODO: 自环等
                pass
            else:
                rplcd = simplify_replace_edges(subgr, sub_edge_points, i, j, path=pd.concat([left_path, right_path]), is_directed=is_directed)
                subgr = rplcd["subgr"]
                sub_edge_points = rplcd["sub_edge_points"]
    subgr.remove_nodes_from(list(keep_v[~keep_v].index))
    return {
        "subgr": subgr,
        "sub_edge_points": sub_edge_points
    }


def simplify_determine_nodes_to_keep(subgr: nx.Graph | nx.DiGraph, is_directed, force_keep=[]):
    # 决定在简化过程中保留哪些节点
    node_list = list(subgr.nodes)
    name_check = pd.Series([i in force_keep for i in node_list], index=node_list)
    loop_check = pd.Series(False, index=node_list)
    loop_check.loc[nx.nodes_with_selfloops(subgr)] = True
    degr_check = []
    if is_directed:
        # 有向图过滤入度1且出度1的节点
        in_degrees = subgr.in_degree
        out_degrees = subgr.out_degree
        for node in node_list:
            in_deg = in_degrees[node]
            out_deg = out_degrees[node]
            if in_deg == 1 and out_deg == 1:
                degr_check.append(False)
            else:
                degr_check.append(True)
    else:
        # 无向图过滤度为2的节点
        degrees = subgr.degree
        for node in node_list:
            deg = degrees[node]
            if deg == 2:
                degr_check.append(False)
            else:
                degr_check.append(True)

    return name_check | loop_check | degr_check  # 或运算符，满足一个条件就保留


def simplify_get_neighbours(subgr, is_directed):
    #  获得图上所有节点的邻居节点序号
    name2id = {name: i for i, name in enumerate(subgr.nodes)}  # name -> id(0,1,2...)
    neighs = {}
    if is_directed:
        neighs["neighs_in"] = []
        neighs["neighs_out"] = []
    else:
        neighs["neighs"] = []
    for node in subgr.nodes:
        if is_directed:
            in_neighbors = [name2id[i] for i in list(subgr.predecessors(node))]  # 前继
            out_neighbors = [name2id[i] for i in list(subgr.successors(node))]  # 后驱
            neighs["neighs_in"] .append(in_neighbors)
            neighs["neighs_out"].append(out_neighbors)
        else:
            neighs["neighs"].append(list(subgr.neighbors(node)))
    return neighs


def simplify_get_i(neighs, v_rem, is_directed):
    # 获取入度节点
    if is_directed:
        return neighs["neighs_in"][v_rem][0]  # 暂时认为只有一个前驱
    else:
        return neighs["neighs"][v_rem][0]


def simplify_get_j(neighs, v_rem, is_directed):
    # 获取出度节点
    if is_directed:
        return neighs["neighs_out"][v_rem][0]  # 暂时认为只有一个后继
    else:
        return neighs["neighs"][v_rem][1]


def simplify_get_next(neighs, v_rem, is_directed, left=False, prev=None):
    # 连续删除节点，prev是上一次删除的节点，v_rem是当前删除的节点
    if is_directed:
        if left:
            # 有向图上不会产生与prev重复的节点
            return neighs["neighs_in"][v_rem][0]  # 入度即向前
        else:
            return neighs["neighs_out"][v_rem][0]  # 出度即向后
    else:
        # 无向图上要除去上一步的节点
        neighs["neighs"].remove(prev)


def simplify_get_edge_points_on_path():
    # TODO: 对于边上点的处理
    pass


def simplify_replace_edges(subgr: nx.Graph | nx.DiGraph, sub_edge_points, i, j, path, is_directed):
    # 添加替换旧边的新边
    node_list = list(subgr.nodes)
    path_len = path["weight"].sum()
    subgr.add_edge(node_list[i], node_list[j], weight=path_len, directed=is_directed)  # 添加新边
    if sub_edge_points is not None:
        # TODO:边上的点处理
        simplify_get_edge_points_on_path()
        sub_edge_points = sub_edge_points
    return {
        "subgr": subgr,
        "sub_edge_points": sub_edge_points
    }


def simplify_get_edge(subgr, i, j):
    node_list = list(subgr.nodes)
    return subgr.edges[(node_list[i], node_list[j])]
