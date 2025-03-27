import pytest
import pandas as pd
import numpy as np
from pydynverse.eval.metric_him import (
    calculate_him,
    get_matched_adjacencies,
    him_distance
)

# 测试样例1：复杂网络（simplify=True）
@pytest.fixture
def complex_networks_simplify_true():
    # 构造 net1（含分支、环及交叉边）
    net1 = pd.DataFrame({
        'from': ['A', 'A', 'B', 'C', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'],
        'to':   ['B', 'C', 'D', 'D', 'E', 'F', 'G', 'F', 'G', 'H', 'D', 'H'],
        'length': [1.0, 1.2, 2.0, 2.0, 1.5, 1.5, 2.5, 1.0, 1.0, 1.8, 2.0, 2.2],
        'directed': [True]*12
    })
    # 构造 net2，与 net1 略有不同（例如增加一条交叉边和不同边长）
    net2 = pd.DataFrame({
        'from': ['A', 'A', 'B', 'C', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'E'],
        'to':   ['B', 'C', 'D', 'D', 'F', 'E', 'G', 'G', 'H', 'I', 'D', 'H', 'I'],
        'length': [1.0, 1.2, 2.1, 2.0, 1.7, 1.4, 2.5, 1.1, 1.0, 1.9, 2.1, 2.2, 1.5],
        'directed': [True]*13
    })
    return net1, net2

def test_complex_topology_simplify_true(complex_networks_simplify_true):
    net1, net2 = complex_networks_simplify_true
    # 此处采用简化过程
    sim = calculate_him(net1, net2, simplify=True, gamma=0.1)
    # 同时获取中间计算的邻接矩阵与 HIM 距离
    adj1, adj2 = get_matched_adjacencies(net1, net2, simplify=True)
    norm_adj1 = adj1 / np.sum(adj1)
    norm_adj2 = adj2 / np.sum(adj2)
    d = him_distance(norm_adj1, norm_adj2, gamma=0.1)
    
    print("Complex topology (simplify=True):")
    print("HIM distance d =", d)
    print("Similarity sim =", sim)
    
    # 期望在简化后两网络仍然有一定差异
    assert d > 0.05
    assert sim < 1.0

# 测试样例2：复杂网络（simplify=False）
@pytest.fixture
def complex_networks_simplify_false():
    # 构造 net1（较复杂网络，保留所有原始细节）
    net1 = pd.DataFrame({
        'from': ['A', 'A', 'B', 'C', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'],
        'to':   ['B', 'C', 'D', 'D', 'E', 'F', 'G', 'F', 'G', 'H', 'D', 'H', 'I'],
        'length': [1.0, 1.2, 2.0, 2.0, 1.5, 1.5, 2.5, 1.0, 1.0, 1.8, 2.0, 2.2, 1.3],
        'directed': [True]*13
    })
    # 构造 net2，与 net1 略有不同：增加多条额外边和权重变化
    net2 = pd.DataFrame({
        'from': ['A', 'A', 'B', 'C', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'E', 'C'],
        'to':   ['B', 'C', 'D', 'D', 'F', 'E', 'G', 'G', 'H', 'I', 'D', 'H', 'I', 'J', 'F'],
        'length': [1.0, 1.2, 2.1, 2.0, 1.7, 1.4, 2.5, 1.1, 1.0, 1.9, 2.1, 2.2, 1.3, 1.8, 1.6],
        'directed': [True]*15
    })
    return net1, net2

def test_complex_topology_simplify_false(complex_networks_simplify_false):
    net1, net2 = complex_networks_simplify_false
    # 保持原始细节，不进行简化
    sim = calculate_him(net1, net2, simplify=False, gamma=0.1)
    # 同时获取中间计算的邻接矩阵与 HIM 距离
    adj1, adj2 = get_matched_adjacencies(net1, net2, simplify=False)
    norm_adj1 = adj1 / np.sum(adj1)
    norm_adj2 = adj2 / np.sum(adj2)
    d = him_distance(norm_adj1, norm_adj2, gamma=0.1)
    
    print("Complex topology (simplify=False):")
    print("HIM distance d =", d)
    print("Similarity sim =", sim)
    
    # 在不简化的情况下，两网络的差异会更明显
    assert d > 0.05
    assert sim < 1.0

if __name__ == "__main__":
    pytest.main(["-v", __file__])
