import networkx as nx
import numpy as np
import pandas as pd
from typing import Tuple, List
from scipy.linalg import eigvalsh

# HIM公式体系完整定义
"""
Hybrid Ipsen-Mikhailov (HIM) 指标体系：

1. HIM_distance = sqrt( (IM² + (γ*Hamming)²) / (1 + γ²) )
2. similarity = max(0, 1 - HIM_distance)

其中各分量计算如下：

█ Ipsen-Mikhailov (IM) 谱距离 █
IM = ||√Λ₁ - √Λ₂||₂ = sqrt(∑(sqrt(λ_i^{(1)}) - sqrt(λ_i^{(2)}))²)

计算步骤：
1. 对邻接矩阵A计算拉普拉斯矩阵：L = D - A，D为度矩阵
2. 计算L的特征值：Λ = {λ_1, λ_2, ..., λ_n}（按升序排列）
3. 对两个网络的Λ₁和Λ₂取平方根后计算欧氏距离

█ Hamming 边权距离 █
Hamming = ||A₁ - A₂||_1 = ∑|A₁[i,j] - A₂[i,j]|

计算步骤：
1. 对两个邻接矩阵对应元素求绝对差
2. 对所有元素的差值求和

█ 参数说明 █
- γ (默认0.1): 调节谱距离与边权距离的权重
- max(0, ...): 确保相似度不小于0
"""
from ..util import random_time_string  
from ..wrap import simplify_networkx_network  

#  边处理函数 
def insert_two_nodes_into_selfloop(df: pd.DataFrame) -> pd.DataFrame:
    """自环边分割处理（严格实现R代码逻辑）"""
    self_loops = df[df['from'] == df['to']].copy()
    if self_loops.empty:
        return df
    
    new_edges = []
    for _, row in self_loops.iterrows():
        new_nodes = [random_time_string() for _ in range(2)]
        split_length = row['length'] / 3
        
        new_edges.extend([
            (row['from'], new_nodes[0], split_length, row['directed']),
            (new_nodes[0], new_nodes[1], split_length, row['directed']),
            (new_nodes[1], row['to'], split_length, row['directed'])
        ])
    
    return pd.concat([
        df[df['from'] != df['to']],
        pd.DataFrame(new_edges, columns=df.columns)
    ], ignore_index=True)

def insert_one_node_into_duplicate_edges(df: pd.DataFrame) -> pd.DataFrame:
    """重复边处理（严格实现R代码逻辑）"""
    edge_keys = df.apply(lambda x: f"{x['from']}#{x['to']}", axis=1)
    dup_mask = edge_keys.isin(edge_keys.value_counts()[edge_keys.value_counts() >= 2].index)
    
    new_edges = []
    for idx in df[dup_mask].index:
        row = df.loc[idx]
        new_node = random_time_string()
        split_length = row['length'] / 2
        
        new_edges.extend([
            (row['from'], new_node, split_length, row['directed']),
            (new_node, row['to'], split_length, row['directed'])
        ])
    
    if new_edges:
        return pd.concat([
            df[~dup_mask],
            pd.DataFrame(new_edges, columns=df.columns)
        ], ignore_index=True)
    return df

def change_single_edge_into_double(df: pd.DataFrame) -> pd.DataFrame:
    """单边分割处理（严格实现R代码逻辑）"""
    if len(df) == 1 and df.iloc[0]['from'] != df.iloc[0]['to']:
        row = df.iloc[0]
        new_node = random_time_string()
        split_length = row['length'] / 2
        
        return pd.DataFrame([
            (row['from'], new_node, split_length, row['directed']),
            (new_node, row['to'], split_length, row['directed'])
        ], columns=df.columns)
    return df

#  网络预处理模块 
def process_simplified_network(
    df: pd.DataFrame,
    directed: bool
) -> pd.DataFrame:
    """
    完整网络处理流程（对应R代码的get_matched_adjacencies中的simplify部分）
    """
    # Step 1: 过滤零长度自环边
    filtered_df = df[(df['from'] != df['to']) | (df['length'] != 0)].copy()
    
    # Step 2: 转换为无向图（强制转换）
    G = nx.from_pandas_edgelist(
        filtered_df.rename(columns={'length': 'weight'}),
        source='from',  # 修复点：指定起点列
        target='to',    # 修复点：指定终点列
        edge_attr=True,
        create_using=nx.Graph  # 强制转换为无向图
    )
    
    # Step 3: 简化网络
    G_simplified = simplify_networkx_network(G)
    
    # Step 4: 转换回DataFrame
    simplified_df = nx.to_pandas_edgelist(G_simplified).rename(columns={
        'source': 'from',  # 强制列名恢复为 from/to
        'target': 'to',
        'weight': 'length'
    })
    simplified_df['directed'] = directed  # 保留原始方向标记
    
    # Step 5: 应用边处理流程
    processed_df = (
        simplified_df.pipe(insert_two_nodes_into_selfloop)
                     .pipe(change_single_edge_into_double)
                     .pipe(insert_one_node_into_duplicate_edges)
    )
    
    return processed_df

#  邻接矩阵处理模块 
def get_adjacency_lengths(df: pd.DataFrame) -> np.ndarray:
    """生成邻接矩阵（严格对齐R代码逻辑）"""
    nodes = sorted(set(df['from']).union(set(df['to'])))
    node_idx = {node: i for i, node in enumerate(nodes)}
    size = len(nodes)
    
    adj = np.zeros((size, size), dtype=np.float64)
    for _, row in df.iterrows():
        i = node_idx[row['from']]
        j = node_idx[row['to']]
        adj[i, j] += row['length']
        if not row['directed']:
            adj[j, i] += row['length']
    return adj

def pad_matrix(mat: np.ndarray, target_size: int) -> np.ndarray:
    """矩阵填充（严格实现R的complete_matrix逻辑）"""
    pad_size = target_size - mat.shape[0]
    if pad_size > 0:
        return np.pad(mat, ((0, pad_size), (0, pad_size)), mode='constant')
    return mat

# 对齐
def get_matched_adjacencies(
    net1: pd.DataFrame,
    net2: pd.DataFrame,
    simplify: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    严格实现R代码的get_matched_adjacencies逻辑
    """
    # 保留原始方向性
    directed1 = net1['directed'].any()
    directed2 = net2['directed'].any()
    
    if simplify:
        processed_net1 = process_simplified_network(net1, directed1)
        processed_net2 = process_simplified_network(net2, directed2)
    else:
        processed_net1 = net1.copy()
        processed_net2 = net2.copy()
    
    # 生成邻接矩阵
    adj1 = get_adjacency_lengths(processed_net1)
    adj2 = get_adjacency_lengths(processed_net2)
    
    # 统一矩阵维度
    max_size = max(adj1.shape[0], adj2.shape[0])
    adj1 = pad_matrix(adj1, max_size)
    adj2 = pad_matrix(adj2, max_size)
    
    return adj1, adj2

#  核心计算模块 
def calculate_him(
    net1: pd.DataFrame,
    net2: pd.DataFrame,
    simplify: bool = True,
    ga: float = 0.1
) -> float:
    """
    严格实现R代码的calculate_him逻辑
    """
    adj1, adj2 = get_matched_adjacencies(net1, net2, simplify)
    
    # 空图检查
    if np.all(adj1 == 0) or np.all(adj2 == 0):
        return 0.0
    
    # 直接计算距离，跳过归一化
    him_distance = compute_him_distance(adj1, adj2, ga)
    return max(0.0, 1.0 - him_distance)

def compute_him_distance(
    m1: np.ndarray,
    m2: np.ndarray,
    ga: float
) -> float:
    """HIM距离核心计算"""
    # 计算Ipsen-Mikhailov距离（修复：特征值排序）
    L1 = compute_laplacian(m1)
    L2 = compute_laplacian(m2)
    eig1 = eigvalsh(L1)  # eigvalsh 默认返回升序
    eig2 = eigvalsh(L2)
    im = np.linalg.norm(np.sqrt(eig1) - np.sqrt(eig2))  # 正确平方根
    
    # 计算Hamming距离
    hamming = np.abs(m1 - m2).sum()
    
    # 组合公式
    return np.sqrt((im**2 + (ga * hamming)**2) / (1 + ga**2))
def compute_laplacian0(matrix: np.ndarray) -> np.ndarray:
    """生成对称拉普拉斯矩阵"""
    matrix_sym = (matrix + matrix.T) / 2  # 强制对称化
    D = np.diag(matrix_sym.sum(axis=1))
    return D - matrix_sym

def compute_laplacian(matrix: np.ndarray) -> np.ndarray:
    """生成拉普拉斯矩阵（无向图）"""
    D = np.diag(matrix.sum(axis=1))
    return D - matrix