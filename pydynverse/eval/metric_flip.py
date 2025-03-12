import numpy as np
import networkx as nx
import pandas as pd
from itertools import combinations
from typing import Union, Tuple, List, Optional
from ..wrap import simplify_networkx_network
from ..util import random_time_string
from math import comb

def insert_two_nodes_into_selfloop(G: nx.Graph) -> nx.Graph:
    """处理自环边：插入两个新节点"""
    edges_to_add = []
    edges_to_remove = []
    
    for u, v, data in G.edges(data=True):
        if u == v:
            edges_to_remove.append((u, v))
            new_node1 = random_time_string()
            new_node2 = random_time_string()
            length = data.get('length', 1.0) / 3
            edges_to_add.extend([
                (u, new_node1, {'length': length}),
                (new_node1, new_node2, {'length': length}),
                (new_node2, u, {'length': length})
            ])
    
    G.remove_edges_from(edges_to_remove)
    G.add_edges_from(edges_to_add)
    return G

def change_single_edge_into_double(G: nx.Graph) -> nx.Graph:
    """处理单一边：插入一个新节点"""
    if len(G.edges()) == 1 and not any(u == v for u, v in G.edges()):
        edges = list(G.edges(data=True))
        G.clear()
        u, v, data = edges[0]
        new_node = random_time_string()
        length = data.get('length', 1.0) / 2
        G.add_edges_from([
            (u, new_node, {'length': length}),
            (new_node, v, {'length': length})
        ])
    return G

def insert_one_node_into_duplicate_edges(G: nx.Graph) -> nx.Graph:
    """处理重复边：插入新节点（支持多重图）"""
    edge_counts = {}
    is_multi = isinstance(G, (nx.MultiGraph, nx.MultiDiGraph))
    is_directed = G.is_directed()

    # 遍历边并统计重复
    if is_multi:
        for u, v, key in G.edges(keys=True):
            edge_key = (u, v) if is_directed else tuple(sorted((u, v)))
            edge_counts[edge_key] = edge_counts.get(edge_key, 0) + 1
    else:
        for u, v in G.edges():
            edge_key = (u, v) if is_directed else tuple(sorted((u, v)))
            edge_counts[edge_key] = edge_counts.get(edge_key, 0) + 1

    edges_to_add = []
    edges_to_remove = []

    for (u, v), count in edge_counts.items():
        if count >= 2:
            # 提取所有原始边数据并计算总长度
            if is_multi:
                edges = []
                for key in G[u][v]:
                    edges.append(G[u][v][key])
            else:
                edges = [G[u][v]] if is_directed else [G[u][v]]
            
            total_length = sum(e.get('length', 1.0) for e in edges)
            avg_length_per_segment = total_length / 2  # 拆分为两条边，总长度保持

            # 删除所有重复边
            if is_multi:
                edges_to_remove.extend([(u, v, key) for key in G[u][v]])
            else:
                edges_to_remove.extend([(u, v)])

            # 插入新边
            new_node = random_time_string()
            edges_to_add.extend([
                (u, new_node, {'length': avg_length_per_segment}),
                (new_node, v, {'length': avg_length_per_segment})
            ])

    # 执行边操作
    G.remove_edges_from(edges_to_remove)
    G.add_edges_from(edges_to_add)
    return G

def _create_adjacency_matrix(G: nx.Graph, nodes: List[str]) -> np.ndarray:
    """创建对称邻接矩阵（严格处理边权）"""
    size = len(nodes)
    adj = np.zeros((size, size), dtype=float)
    node_index = {n: i for i, n in enumerate(nodes)}
    
    for u, v, data in G.edges(data=True):
        i = node_index[u]
        j = node_index[v]
        adj[i][j] += data.get('length', 1.0)
        adj[j][i] = adj[i][j]  # 保证对称性
        
    return adj

def get_matched_adjacencies(
    net1: nx.Graph,
    net2: nx.Graph,
    simplify: bool = False
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """生成无向图邻接矩阵并对齐"""
    def process_network(net: nx.Graph) -> nx.Graph:
        """必须执行的处理步骤（无论是否简化）"""
        net = insert_two_nodes_into_selfloop(net)
        net = change_single_edge_into_double(net)
        net = insert_one_node_into_duplicate_edges(net)
        return net

    # Step 1: 网络简化
    if simplify:
        net1 = simplify_networkx_network(net1)
        net2 = simplify_networkx_network(net2)

    # Step 2: 必须执行的处理流程(处理自环，重边，单一边)
    net1 = process_network(net1)
    net2 = process_network(net2)

    # 合并并排序所有节点
    #all_nodes = sorted(set(net1.nodes()) | set(net2.nodes()))#(上一版本)
    # 合并并排序所有节点时保留原始节点类型
    all_nodes = sorted(set(net1.nodes()) | set(net2.nodes()), 
                      key=lambda x: (isinstance(x, str), str(x)))
    # 创建邻接矩阵
    adj1 = _create_adjacency_matrix(net1, all_nodes)
    adj2 = _create_adjacency_matrix(net2, all_nodes)

    # 矩阵维度对齐（补充空节点）
    """max_size = max(adj1.shape[0], adj2.shape[0])
    adj1 = np.pad(adj1, ((0, max_size-adj1.shape[0]), (0, max_size-adj1.shape[0])))
    adj2 = np.pad(adj2, ((0, max_size-adj2.shape[0]), (0, max_size-adj2.shape[0])))
"""
    max_size = len(all_nodes)  # 直接使用合并后的节点数作为统一维度
    adj1.resize((max_size, max_size), refcheck=False)
    adj2.resize((max_size, max_size), refcheck=False)
    return adj1, adj2, all_nodes

def calculate_edge_membership0(adj: np.ndarray) -> np.ndarray:
    """精确实现R的边-节点关系矩阵（过滤0权边）"""
    #草案，可能废除

    n = adj.shape[0]
    tril_indices = np.tril_indices(n, k=-1)
    edge_mask = adj[tril_indices] > 0  #边点映射只筛选出有边的地方
    
    # 仅处理实际存在的边
    valid_edges = np.where(edge_mask)[0]
    edge_count = len(valid_edges)
    
    membership = np.zeros((edge_count, n), dtype=int)
    for valid_idx, idx in enumerate(valid_edges):  # 关键修改：遍历有效边的索引
        i = tril_indices[0][idx]
        j = tril_indices[1][idx]
        membership[valid_idx, i] = 1
        membership[valid_idx, j] = 1
    return membership

def calculate_edge_membership(adj: np.ndarray) -> np.ndarray:
    """生成包含所有可能边的成员矩阵
    
    这里生成矩阵的维度是(n_edges,n),其实n就相当于是网络邻接矩阵中蕴含的节点数
    含义：边的编号是[0,n_edges),其实可以对应矩阵的每一行，然后对于每一行，如果
    行内全为0，则说明对应的边没有任何连接的节点信息，也就是说不存在这条边。
    请注意，n_edges是矩阵下三角的元素数量，也就是说这个时候它其实就相当于一个
    无向图可能存在的边的数量。(但是这里需要考虑一些问题，我是否需要把0权边加上？或者说把空边也纳入到矩阵里？)
    """
    n = adj.shape[0]
    tril_indices = np.tril_indices(n, k=-1)
    n_edges = len(tril_indices[0])
    membership = np.zeros((n_edges, n), dtype=int)
    for idx in range(n_edges):
        i, j = tril_indices[0][idx], tril_indices[1][idx]
        # 无论原来是否存在边，都记录该边连接的节点信息
        membership[idx, i] = 1
        membership[idx, j] = 1
    return membership

def check_degrees_max(degree_vectors: np.ndarray, target_max: int) -> np.ndarray:
    """最大度检查（向量化实现）"""
    return np.max(degree_vectors, axis=1) == target_max

def check_degrees_min(degree_vectors: np.ndarray, target_min: int) -> np.ndarray:
    """最小度检查（向量化实现）"""
    return np.min(degree_vectors, axis=1) == target_min

def check_degrees_sorted(degree_vectors: np.ndarray, target_sorted: np.ndarray) -> np.ndarray:
    """排序度检查（精确匹配）"""
    return np.all(np.sort(degree_vectors, axis=1) == target_sorted, axis=1)

def generate_edge_flip_vectors(flips: np.ndarray, adj: np.ndarray, tril_indices: tuple) -> np.ndarray:
    """完整处理各种边界情况的翻转向量生成"""
    edge_vector = adj[tril_indices].copy()
    
    # 处理空输入
    if flips.size == 0:
        return np.empty((0, edge_vector.size), dtype=int)
    
    # 统一处理不同维度输入
    if flips.dtype == object:
        # 处理不规则长度的组合
        flips_list = [f.astype(int) for f in flips.ravel() if f.size > 0]
        if not flips_list:
            return np.empty((0, edge_vector.size), dtype=int)
        
        # 转换所有元素为整数
        max_len = max(len(f) for f in flips_list)
        padded_flips = np.full((len(flips_list), max_len), -1, dtype=int)
        for i, f in enumerate(flips_list):
            padded_flips[i, :len(f)] = f
        
        # 过滤无效索引
        valid_mask = (padded_flips >= 0) & (padded_flips < edge_vector.size)
        flips = np.where(valid_mask, padded_flips, -1)
    else:
        # 转换标准数组
        flips = flips.astype(int)
    
    # 处理无效索引(将数组元素限制在[0,edge_vector.size-1])
    flips = np.clip(flips, 0, edge_vector.size-1)
    
    # 生成向量
    vectors = np.tile(edge_vector, (flips.shape[0], 1))
    np.put_along_axis(vectors, flips, 1 - np.take_along_axis(vectors, flips, axis=1), axis=1)
    return vectors

def process_batch(batch, adj1, adj2, tril_indices, edge_membership):
    """完整处理各种边界情况的批次处理"""
    if not batch:
        return False
    
    # 转换批次数据
    try:
        # 尝试转换为标准整数数组
        batch_array = np.stack(batch).astype(int)
    except (ValueError, TypeError):
        # 处理不规则长度的组合
        batch_array = np.array([np.array(f, dtype=int) for f in batch], dtype=object)
    
    # 处理全空批次
    if batch_array.size == 0:
        return False
    
    # 生成翻转向量
    try:
        flip_vectors = generate_edge_flip_vectors(batch_array, adj1, tril_indices)
    except Exception as e:
        print(f"Error generating flip vectors: {e}")
        return False
    
    # 计算度数向量
    degree_vectors = flip_vectors @ edge_membership
    
    # 多级过滤
    candidates = np.arange(len(batch))
    
    # 1. 最大度过滤
    target_max = np.max(np.sum(adj2, axis=1))
    candidates = candidates[check_degrees_max(degree_vectors, target_max)]
    
    # 2. 最小度过滤
    if candidates.size > 0:
        target_min = np.min(np.sum(adj2, axis=1))
        candidates = candidates[check_degrees_min(degree_vectors[candidates], target_min)]
    
    # 3. 排序度过滤
    if candidates.size > 0:
        sorted_target = np.sort(np.sum(adj2, axis=1))
        candidates = candidates[check_degrees_sorted(degree_vectors[candidates], sorted_target)]
    
    # 4. 同构检查
    if candidates.size > 0:
        target_graph = nx.from_numpy_array(adj2)
        # Remove isolated nodes from target_graph
        target_graph.remove_nodes_from(list(nx.isolates(target_graph)))
        for idx in candidates:
            try:
                flip = batch[idx]
                new_adj = flip_adj(flip, adj1.copy(), tril_indices)
                candidate_graph = nx.from_numpy_array(new_adj)
                # Remove isolated nodes from candidate_graph
                candidate_graph.remove_nodes_from(list(nx.isolates(candidate_graph)))
                if nx.is_isomorphic(candidate_graph, target_graph):
                    global newadj1
                    newadj1 = new_adj
                    return True
            except Exception as e:
                print(f"Error in isomorphism check: {e}")
                continue
    return False

def flip_adj0(flip_indices: np.ndarray, adj: np.ndarray, tril_indices: tuple) -> np.ndarray:
    """仅翻转 candidate 指定位置的邻接矩阵，并保持对称性"""
    new_adj = adj.copy()
    lower_vals = new_adj[tril_indices].copy()
    # 显式对每个候选索引进行翻转
    for idx in flip_indices:
        lower_vals[idx] = 1 - lower_vals[idx]
    new_adj[tril_indices] = lower_vals
    return np.maximum(new_adj, new_adj.T)

def flip_adj(flip_indices: np.ndarray, adj: np.ndarray, tril_indices: tuple) -> np.ndarray:
    new_adj = adj.copy()
    lower_vals = new_adj[tril_indices].copy()
    # 扁平化候选索引，确保得到的是一维数组
    for idx in np.ravel(flip_indices):
        lower_vals[idx] = 1 - lower_vals[idx]
    new_adj[tril_indices] = lower_vals
    return np.maximum(new_adj, new_adj.T)

def combn_nice(elements, k: int) -> list:
    """稳健的组合生成（兼容空集）"""
    if k <= 0 or not elements:
        return [np.array([], dtype=int)]
    return [np.array(c) for c in combinations(elements, k)]

def calculate_edge_flip(
    net1: nx.Graph,
    net2: nx.Graph,
    return_type: str = "score",
    simplify: bool = False,
    limit_flips: int = 5,
    limit_combinations: int = comb(25, 4)
) -> Union[float, dict]:
    """完整实现边翻转分数计算"""
    #如果图同构，那直接返回1
    if nx.is_isomorphic(net1, net2):
        return 1
    
    # 对齐邻接矩阵
    adj1, adj2, node_order = get_matched_adjacencies(net1, net2, simplify)
    
    # 转换为二进制邻接矩阵
    adj1_bin = (adj1 > 0).astype(int)
    adj2_bin = (adj2 > 0).astype(int)
    
    # 下三角处理（严格对应R的lower.tri）
    n = adj1_bin.shape[0]
    tril_indices = np.tril_indices(n, k=-1)
    
    # 计算边差异和两图边总数。
    # 使用带符号的差值
    diff = int(np.sum(adj2_bin[tril_indices]) - np.sum(adj1_bin[tril_indices]))
    edge_diff = abs(diff)
    upper_bound = max(np.sum(adj1_bin[tril_indices]) + np.sum(adj2_bin[tril_indices]), 1)


    #搜索参数和结果参数的初始化。
    found = False
    newadj1 = None
    #获得一个边点映射矩阵，每一行代表可能的边，行中两个有值的点的位置对(i,j)代表边连接的节点编号。
    edge_membership = calculate_edge_membership(adj1_bin)
    
    # 初始化翻转次数
    n_flips = max(abs(edge_diff) - 2, 0)
    if edge_diff % 2 != 0:
        n_flips += 1

    while True:
        current_flips = n_flips
        if current_flips > min(upper_bound, limit_flips):
            break
            # 根据 diff 的正负决定 n_add 与 n_remove
        if diff >= 0:
            n_add = (current_flips + diff) // 2
            n_remove = current_flips - n_add
        else:
            n_remove = (current_flips - diff) // 2
            n_add = current_flips - n_remove
        
        # 有效性检查
        if any([n_add < 0, n_remove < 0, (n_add + n_remove) != current_flips]):
            n_flips += 2
            continue
        
        # 生成候选边操作
        possible_add = np.where(adj1_bin[tril_indices] == 0)[0].tolist()
        possible_remove = np.where(adj1_bin[tril_indices] == 1)[0].tolist()
        
        # 生成可能的组合数
        add_comb = comb(len(possible_add), n_add) if n_add > 0 else 1
        remove_comb = comb(len(possible_remove), n_remove) if n_remove > 0 else 1
        if (add_comb > limit_combinations) or (remove_comb > limit_combinations):
            n_flips += 2
            continue
    
        # 生成组合（优化内存使用）
        add_combs = combn_nice(possible_add, n_add)
        remove_combs = combn_nice(possible_remove, n_remove)
        
        # 新增：检查笛卡尔积总数
        if len(add_combs) * len(remove_combs) > limit_combinations:
            n_flips += 2
            continue
        
        # 生成笛卡尔积（使用生成器优化）
        edge_flips = (
            np.concatenate([add, remove])
            for add in add_combs
            for remove in remove_combs
        )
        
        # 分块处理（提升性能,考虑是否需要循环？）
        batch_size = 1000
        current_batch = []
        
        for flip in edge_flips:
            current_batch.append(flip)
            if len(current_batch) == batch_size:
                if process_batch(current_batch, adj1_bin, adj2_bin, 
                               tril_indices, edge_membership):
                    found = True
                    break
                current_batch = []
        
        # 处理剩余批次
        if not found and current_batch:
            if process_batch(current_batch, adj1_bin, adj2_bin,
                           tril_indices, edge_membership):
                found = True
        
        if found:
            break
        else:
            n_flips += 2
    
    if not found:
        n_flips = upper_bound
    
    score = 1 - n_flips / upper_bound
    return {
        "score": score,
        "newadj1": newadj1,
        "oldadj1": adj1_bin,
        "node_order": node_order
    } if return_type == "all" else score

