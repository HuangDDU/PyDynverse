import numpy as np
import networkx as nx
import pandas as pd
from itertools import combinations, product
from typing import Union, List, Tuple, Dict
from math import comb
from ..wrap import simplify_networkx_network
from ..util import random_time_string


def calculate_edge_flip(
    net1: nx.Graph,
    net2: nx.Graph,
    return_type: str = "score",
    simplify: bool = False,
    limit_flips: int = 5,
    limit_combinations: int = 12650
) -> Union[float, dict]:
    # 获取对齐后的邻接矩阵
    adj1, adj2 = get_matched_adjacencies(net1, net2, simplify)
    
    # 转换为二进制邻接矩阵
    adj1_bin = (adj1 > 0).astype(int)
    adj2_bin = (adj2 > 0).astype(int)
    
    # 计算边差异
    triu_indices = np.triu_indices_from(adj1_bin, k=1)
    edge_diff = adj2_bin[triu_indices].sum() - adj1_bin[triu_indices].sum()
    
    # 计算upper bound
    upper_bound = adj1_bin[triu_indices].sum() + adj2_bin[triu_indices].sum()
    upper_bound = max(upper_bound, 1) if upper_bound > 0 else 1

    # 初始化搜索参数
    found = False
    n_flips = abs(edge_diff) - 2 if edge_diff != 0 else 0
    newadj1 = None

    # 准备边成员关系矩阵
    edge_membership = calculate_edge_membership(adj1_bin)

    # 主搜索循环
    while not found and n_flips <= upper_bound:
        n_flips += 2
        
        if n_flips > limit_flips:
            n_flips = upper_bound
            break
            
        # 计算需要添加/删除的边数
        n_add = (n_flips + edge_diff) // 2
        n_remove = (n_flips - edge_diff) // 2
        
        # 生成可能的边操作组合
        possible_add = np.where(adj1_bin[triu_indices] == 0)[0].tolist()
        possible_remove = np.where(adj1_bin[triu_indices] > 0)[0].tolist()
        
        # 有效性检查
        if n_add < 0 or n_remove < 0:
            continue
            
        # 生成组合
        add_combs = combn_nice(possible_add, n_add)
        remove_combs = combn_nice(possible_remove, n_remove)
        
        # 组合矩阵生成
        if n_add > 0 and n_remove > 0:
            edge_flips = np.vstack([
                np.repeat(add_combs, remove_combs.shape[1], axis=1),
                np.tile(remove_combs, (1, add_combs.shape[1]))
            ])
        elif n_add > 0:
            edge_flips = add_combs
        else:
            edge_flips = remove_combs

        # 分块处理（每1000个组合）
        for chunk in np.array_split(edge_flips, max(1, edge_flips.shape[1]//1000), axis=1):
            # 生成邻接矩阵向量
            flip_vectors = generate_edge_flip_vectors(chunk, adj1_bin)
            
            # 计算度数向量
            degree_vectors = flip_vectors.T @ edge_membership
            
            # 四层过滤
            selected = np.arange(degree_vectors.shape[0])
            
            # 1. 最大度过滤
            max_check = degree_vectors.max(axis=1) == adj2_bin.sum(axis=1).max()
            selected = selected[max_check]
            if selected.size == 0:
                continue
                
            # 2. 最小度过滤
            min_check = degree_vectors[selected].min(axis=1) == adj2_bin.sum(axis=1).min()
            selected = selected[min_check]
            if selected.size == 0:
                continue
                
            # 3. 排序度过滤
            sorted_check = np.all(
                np.sort(degree_vectors[selected], axis=1) == np.sort(adj2_bin.sum(axis=1)), 
                axis=1
            )
            selected = selected[sorted_check]
            if selected.size == 0:
                continue
                
            # 4. 同构检查
            for idx in selected:
                new_adj = flip_adj(chunk[:, idx], adj1_bin.copy())
                if nx.is_isomorphic(
                    nx.from_numpy_array(new_adj), 
                    nx.from_numpy_array(adj2_bin),
                    edge_match=lambda e1,e2: e1['weight']==e2['weight']
                ):
                    found = True
                    newadj1 = new_adj
                    break
            if found:
                break

    # 结果处理
    if not found:
        raise RuntimeError("No valid mapping found")
    
    score = 1 - n_flips / upper_bound
    
    if return_type == "all":
        return {
            "score": score,
            "newadj1": newadj1,
            "oldadj1": adj1_bin
        }
    else:
        return score


def get_matched_adjacencies(net1: nx.Graph, net2: nx.Graph, simplify: bool) -> Tuple[np.ndarray, np.ndarray]:
    """对齐邻接矩阵并处理特殊边"""
    if simplify:
        net1 = simplify_networkx_network(net1)
        net2 = simplify_networkx_network(net2)
    
    # 获取所有节点并排序
    all_nodes = sorted(set(net1.nodes()) | set(net2.nodes()))
    
    # 确保节点存在
    for node in all_nodes:
        if node not in net1.nodes():
            net1.add_node(node)
        if node not in net2.nodes():
            net2.add_node(node)
    
    # 生成邻接矩阵
    adj1 = nx.to_numpy_array(net1, nodelist=all_nodes, weight='length', nonedge=0)
    adj2 = nx.to_numpy_array(net2, nodelist=all_nodes, weight='length', nonedge=0)
    
    return adj1, adj2


def process_special_edges(net: nx.Graph, adj: np.ndarray, nodes: list) -> np.ndarray:
    """处理自环边和重复边"""
    # 插入自环边处理
    for u, v in net.edges():
        if u == v:
            new_nodes = [random_time_string() for _ in range(2)]
            adj = insert_nodes_into_edge(adj, nodes, u, v, new_nodes)
    
    # 处理重复边
    edge_counts = pd.Series(net.edges()).value_counts()
    for (u, v), cnt in edge_counts.items():
        if cnt > 1:
            new_node = random_time_string()
            adj = insert_nodes_into_edge(adj, nodes, u, v, [new_node])
    
    return adj

def insert_nodes_into_edge(adj: np.ndarray, nodes: list, u: str, v: str, new_nodes: list) -> np.ndarray:
    """插入节点到边中"""
    u_idx = nodes.index(u)
    v_idx = nodes.index(v)
    
    # 扩展邻接矩阵
    for n in new_nodes:
        if n not in nodes:
            nodes.append(n)
            adj = np.pad(adj, [(0,1), (0,1)])
    
    # 重建边关系
    prev_node = u
    for n in new_nodes:
        n_idx = nodes.index(n)
        adj[u_idx, n_idx] = adj[n_idx, u_idx] = adj[u_idx, v_idx]/len(new_nodes)+1
        prev_node = n
    
    adj[prev_node, v_idx] = adj[v_idx, prev_node] = adj[u_idx, v_idx]/len(new_nodes)+1
    adj[u_idx, v_idx] = adj[v_idx, u_idx] = 0
    
    return adj

def calculate_edge_membership(adj: np.ndarray) -> np.ndarray:
    """计算边-节点关系矩阵（修复shape mismatch问题）"""
    triu_indices = np.triu_indices_from(adj, k=1)
    triu_size = len(triu_indices[0])  # 直接使用索引数量
    
    edge_mapper = np.zeros_like(adj, dtype=int)
    edge_mapper[triu_indices] = np.arange(triu_size)
    
    membership = []
    for i in range(adj.shape[0]):
        row = ((edge_mapper[i,:] >= 0) | (edge_mapper[:,i] >= 0)).astype(int)
        membership.append(row)
    
    return np.column_stack(membership)

def generate_edge_flip_vectors(edge_flips: np.ndarray, adj: np.ndarray) -> np.ndarray:
    """生成边翻转向量"""
    triu_flat = adj[np.triu_indices_from(adj, k=1)].copy()
    for flip in edge_flips.T:
        triu_flat[flip] = 1 - triu_flat[flip]
    return triu_flat.reshape(1, -1)

def flip_adj(flip_indices: np.ndarray, adj: np.ndarray) -> np.ndarray:
    """应用翻转得到新邻接矩阵"""
    new_adj = adj.copy()
    triu_indices = np.triu_indices_from(adj, k=1)
    new_adj[triu_indices][flip_indices] = 1 - new_adj[triu_indices][flip_indices]
    return np.maximum(new_adj, new_adj.T)

def combn_nice(elements, k):
    """Generate all combinations of k elements from the list, returns a numpy array."""
    if k == 0:
        return np.empty((0, 0), dtype=int)
    if len(elements) < k:
        return np.empty((0, k), dtype=int)
    combs = list(combinations(elements, k))
    return np.array(combs, dtype=int) if combs else np.empty((0, k), dtype=int)
