import pytest
import pandas as pd
import numpy as np
from pydynverse.eval.metric_him import calculate_him, get_matched_adjacencies, him_distance


# 构造相似的树状网络数据
@pytest.fixture
def similar_trees():
    # 构造树形网络 net1
    # 树结构：根节点 "A"，下面三个分支 "B", "C", "D"；每个分支再扩展两个子节点
    net1_edges = [
        ('A', 'B', 1.0), ('A', 'C', 1.0), ('A', 'D', 1.0),
        ('B', 'B1', 0.8), ('B', 'B2', 0.9),
        ('C', 'C1', 0.7), ('C', 'C2', 0.85),
        ('D', 'D1', 0.95), ('D', 'D2', 0.8)
    ]
    net1 = pd.DataFrame(net1_edges, columns=['from', 'to', 'length'])
    net1['directed'] = True

    # 构造相似的树形网络 net2
    # 与 net1 基本相同，仅在少数边上略有差异，例如边长度稍有调整，或者增加一条额外边
    net2_edges = [
        ('A', 'B', 1.0), ('A', 'C', 1.0), ('A', 'D', 1.0),
        ('B', 'B1', 0.8), ('B', 'B2', 1.0),   # B2 边长度从0.9变为1.0
        ('C', 'C1', 0.7), ('C', 'C2', 0.85),
        ('D', 'D1', 0.95), ('D', 'D2', 0.8),
        # 增加一条额外的边：C->B，增加轻微的交叉关系
        ('C', 'B', 0.5)
    ]
    net2 = pd.DataFrame(net2_edges, columns=['from', 'to', 'length'])
    net2['directed'] = True

    return net1, net2


# 构造差异较大的树状网络数据
@pytest.fixture
def different_trees():
    # 构造树形网络 net1
    # 使用与上面类似的树状结构
    net1_edges = [
        ('A', 'B', 1.0), ('A', 'C', 1.0), ('A', 'D', 1.0),
        ('B', 'B1', 0.8), ('B', 'B2', 0.9),
        ('C', 'C1', 0.7), ('C', 'C2', 0.85),
        ('D', 'D1', 0.95), ('D', 'D2', 0.8)
    ]
    net1 = pd.DataFrame(net1_edges, columns=['from', 'to', 'length'])
    net1['directed'] = True

    # 构造差异较大的树形网络 net2
    # 这里修改网络结构：改变分支连接和边权，令 net2 与 net1 有较大差异
    net2_edges = [
        ('A', 'B', 1.0), ('A', 'E', 1.2),  # 不再直接连接 A->C，而是 A->E
        ('B', 'B1', 1.5),  ('B', 'B2', 1.4),
        ('E', 'C', 0.9),   ('E', 'F', 1.1),  # E 分支出两个子节点，分别连接 C 和 F
        ('C', 'C1', 0.7), ('C', 'C2', 0.85),
        ('F', 'D', 1.3),  ('F', 'D1', 1.2)   # F 分支再连接 D 和 D1，而非直接从 A->D
    ]
    net2 = pd.DataFrame(net2_edges, columns=['from', 'to', 'length'])
    net2['directed'] = True

    return net1, net2


# 测试：相似树状网络（简化过程可选，这里我们用 simplify=True 进行测试）
def test_similar_trees(similar_trees):
    net1, net2 = similar_trees
    # 使用简化过程（如果你希望观察简化后的结果）
    sim = calculate_him(net1, net2, simplify=True, gamma=0.1)
    # 同时提取中间计算的邻接矩阵和 HIM 距离
    adj1, adj2 = get_matched_adjacencies(net1, net2, simplify=True)
    norm_adj1 = adj1 / np.sum(adj1)
    norm_adj2 = adj2 / np.sum(adj2)
    d = him_distance(norm_adj1, norm_adj2, gamma=0.1)
    
    print("Similar trees (simplify=True):")
    print("HIM distance d =", d)
    print("Similarity sim =", sim)
    
    # 预期相似度较高，相应 HIM 距离应较小（例如小于0.1，根据数据具体结果可能略有波动）
    assert sim > 0.9
    assert d < 0.1

# 测试：差异较大的树状网络（这里用 simplify=False 保留更多原始细节）
def test_different_trees(different_trees):
    net1, net2 = different_trees
    sim = calculate_him(net1, net2, simplify=False, gamma=0.1)
    adj1, adj2 = get_matched_adjacencies(net1, net2, simplify=False)
    norm_adj1 = adj1 / np.sum(adj1)
    norm_adj2 = adj2 / np.sum(adj2)
    d = him_distance(norm_adj1, norm_adj2, gamma=0.1)
    
    print("Different trees (simplify=False):")
    print("HIM distance d =", d)
    print("Similarity sim =", sim)
    
    # 预期差异较大，相似度较低，HIM 距离较大（例如 HIM 距离 > 0.2，相似度 < 0.8）
    assert d > 0.2
    assert sim < 0.8

if __name__ == "__main__":
    pytest.main(["-v", __file__])
