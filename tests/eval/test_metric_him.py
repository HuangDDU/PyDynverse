import pytest
import pandas as pd
import numpy as np
import networkx as nx
from pydynverse.wrap.simplify_networkx_network import simplify_networkx_network 
from pydynverse.util.random_time_string import random_time_string

from pydynverse.eval.metric_him import (
    insert_two_nodes_into_selfloop,
    insert_one_node_into_duplicate_edges,
    change_single_edge_into_double,
    get_adjacency_lengths,
    complete_matrix,
    get_matched_adjacencies,
    laplacian_matrix,
    ipsen_mikhailov_distance_eigen,
    hamming_distance,
    him_distance,
    calculate_him
)

# 构造简单测试网络数据供后续使用
@pytest.fixture
def sample_networks():
    net1 = pd.DataFrame({
        'from': ['M1', 'M2', 'M3'],
        'to':   ['M2', 'M3', 'M3'],  # M3->M3 为自环
        'length': [1.0, 2.0, 0.0],
        'directed': [True, True, True]
    })
    net2 = pd.DataFrame({
        'from': ['M1', 'M2', 'M3'],
        'to':   ['M2', 'M1', 'M1'],
        'length': [1.0, 1.5, 1.0],
        'directed': [True, True, True]
    })
    return net1, net2

@pytest.fixture
def sample_edges_df():
    # 用于测试辅助函数的边表数据
    return pd.DataFrame({
        'from': ['A', 'B', 'C', 'C'],
        'to':   ['A', 'C', 'D', 'D'],
        'length': [3.0, 4.0, 5.0, 5.0],
        'directed': [False, False, False, False]
    })

def test_insert_two_nodes_into_selfloop():
    df = pd.DataFrame({
        'from': ['X', 'Y'],
        'to': ['X', 'Z'],
        'length': [9.0, 2.0],
        'directed': [False, False]
    })
    df_result = insert_two_nodes_into_selfloop(df)
    # 原来自环边 ('X'->'X') 被拆分成两条边，行数应增加
    assert len(df_result) > len(df)
    # 检查拆分后结果中不再存在 from == to 的行
    assert not any(df_result['from'] == df_result['to'])

def test_insert_one_node_into_duplicate_edges(sample_edges_df):
    df_result = insert_one_node_into_duplicate_edges(sample_edges_df)
    # 原始 df 中 'C'->'D' 出现两次，拆分后应新增边，行数增加
    assert len(df_result) > len(sample_edges_df)
    # 检查是否产生了不在原始节点集合内的新节点（随机生成）
    orig_nodes = set(sample_edges_df['from']).union(set(sample_edges_df['to']))
    new_nodes = set(df_result['from']).union(set(df_result['to']))
    assert len(new_nodes - orig_nodes) >= 1

def test_change_single_edge_into_double():
    df_single = pd.DataFrame({
        'from': ['P'],
        'to': ['Q'],
        'length': [6.0],
        'directed': [True]
    })
    df_result = change_single_edge_into_double(df_single)
    # 根据示例逻辑，拆分后应有两条边
    assert len(df_result) == 2

def test_get_adjacency_lengths(sample_edges_df):
    A = get_adjacency_lengths(sample_edges_df)
    # 节点集合为 A, B, C, D，矩阵尺寸应为 4x4
    assert A.shape == (4, 4)
    # 查找 'C' 到 'D' 的值是否为 5.0（注意重复边时后面可能覆盖前面，此处仅测试基本功能）
    nodes = sorted(set(sample_edges_df['from']).union(sample_edges_df['to']))
    idx_C = nodes.index('C')
    idx_D = nodes.index('D')
    np.testing.assert_almost_equal(A[idx_C, idx_D], 5.0)

def test_complete_matrix():
    mat = np.array([[1, 2], [3, 4]])
    new_mat = complete_matrix(mat, 3, fill=0)
    assert new_mat.shape == (3, 3)
    # 检查右下角应为填充值 0
    assert new_mat[2, 2] == 0

def test_laplacian_matrix():
    A = np.array([[0, 1], [1, 0]])
    L = laplacian_matrix(A)
    expected = np.array([[1, -1], [-1, 1]])
    np.testing.assert_allclose(L, expected)

def test_ipsen_mikhailov_distance_eigen():
    A1 = np.array([[0, 1], [1, 0]])
    A2 = np.array([[0, 0.5], [0.5, 0]])
    d = ipsen_mikhailov_distance_eigen(A1, A2)
    # 两个图较接近，距离应大于0且较小
    assert d > 0
    assert d < 1

def test_hamming_distance():
    A1 = np.array([[0, 1], [1, 0]])
    A2 = np.array([[0, 0.5], [0.5, 0]])
    h = hamming_distance(A1, A2)
    # 差异为0.5处各出现两次，总和应为1
    np.testing.assert_almost_equal(h, 1.0)

def test_him_distance():
    # 使用 3x3 矩阵，使归一化后依然不同
    A1 = np.array([
        [0, 1, 1],
        [1, 0, 0],
        [1, 0, 0]
    ])
    A2 = np.array([
        [0, 0.5, 0.2],
        [0.5, 0,   0],
        [0.2, 0,   0]
    ])
    norm_A1 = A1 / A1.sum()
    norm_A2 = A2 / A2.sum()
    d = him_distance(norm_A1, norm_A2, gamma=0.1)
    # HIM 距离应大于0
    assert d > 0

def test_calculate_him(sample_networks):
    net1, net2 = sample_networks
    sim = calculate_him(net1, net2, simplify=True, gamma=0.1)
    # 返回值应在 0 到 1 之间
    assert 0 <= sim <= 1
    # 便于调试时输出结果
    print("HIM similarity:", sim)

# 复杂网络测试：构造两个拓扑差异明显的网络
@pytest.fixture
def complex_networks():
    # net1: 链式结构 A->B->C->D->E
    net1 = pd.DataFrame({
        'from': ['A', 'B', 'C', 'D'],
        'to':   ['B', 'C', 'D', 'E'],
        'length': [1.0, 2.0, 1.0, 2.0],
        'directed': [True, True, True, True]
    })
    
    # net2: 环结构，并有额外边 D->A
    net2 = pd.DataFrame({
        'from': ['A', 'B', 'C', 'D', 'D'],
        'to':   ['B', 'C', 'D', 'E', 'A'],
        'length': [1.0, 1.5, 1.0, 2.0, 3.0],
        'directed': [True, True, True, True, True]
    })
    return net1, net2

def test_complex_topology(complex_networks):
    net1, net2 = complex_networks
    sim = calculate_him(net1, net2, simplify=False, gamma=0.1)
    # 计算 HIM 距离 d 也可以通过中间步骤来查看
    adj1, adj2 = get_matched_adjacencies(net1, net2, simplify=False)
    # 归一化邻接矩阵
    norm_adj1 = adj1 / np.sum(adj1)
    norm_adj2 = adj2 / np.sum(adj2)
    d = him_distance(norm_adj1, norm_adj2, gamma=0.1)
    
    print("Complex topology test:")
    print("HIM distance d =", d)
    print("Similarity sim =", sim)
    
    # 期望 d 应该明显大于 0，不至于太小（具体数值依赖于简化接口和网络归一化效果）
    assert d > 0.05
    # 相似度 sim 则会小于 1，且反映一定差异
    assert sim < 1.0

if __name__ == "__main__":
    pytest.main(["-v", __file__])
