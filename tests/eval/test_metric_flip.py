# test_metric_flip.py
import pytest
import pandas as pd
import networkx as nx
from itertools import product
from typing import Dict, Tuple
import pydynverse as pdv
from pydynverse.eval.metric_flip import *


#函数功能测试
def test_insert_one_node_into_duplicate_edges():
    """验证重复边被合并为两条新边，长度取平均"""
    # 构造含重复边的测试数据（使用MultiDiGraph）
    base_net = pd.DataFrame([
        ["a", "b", 2.0, True],
        ["a", "b", 4.0, True]  # 重复边
    ], columns=["from", "to", "length", "directed"])
    
    # 创建多重有向图
    net = nx.from_pandas_edgelist(base_net, "from", "to", edge_attr=True, create_using=nx.MultiDiGraph())
    processed_net = insert_one_node_into_duplicate_edges(net)
    
    # 断言检查
    assert len(processed_net.edges()) == 2, "应合并为2条新边"
    assert len(processed_net.nodes()) == 3, "应新增1个中间节点"
    edges = list(processed_net.edges(data=True))
    assert all(data["length"] == 3.0 for _, _, data in edges), "每条边长度应为平均值3.0"

def test_change_single_edge_into_double():
    """验证单一边被拆分为两条边，中间插入新节点"""
    # 构造单一边测试数据
    base_net = pd.DataFrame([
        ["a", "b", 2.0, False]
    ], columns=["from", "to", "length", "directed"])
    
    net = nx.from_pandas_edgelist(base_net, "from", "to", edge_attr=True, create_using=nx.Graph())
    processed_net = change_single_edge_into_double(net)
    
    # 断言检查
    assert len(processed_net.edges()) == 2, "应拆分为2条边"
    assert len(processed_net.nodes()) == 3, "应新增1个中间节点"
    edges = list(processed_net.edges(data=True))
    assert edges[0][2]["length"] == 1.0 and edges[1][2]["length"] == 1.0, "每条边长度应为2.0/2=1.0"

def test_insert_two_nodes_into_selfloop():
    """验证自环边被替换成3条新边，且长度正确分割"""
    # 构造含自环边的测试数据
    base_net = pd.DataFrame([
        ["a", "a", 3.0, False]  # 自环边
    ], columns=["from", "to", "length", "directed"])
    
    net = nx.from_pandas_edgelist(base_net, "from", "to", edge_attr=True, create_using=nx.Graph())
    processed_net = insert_two_nodes_into_selfloop(net)
    
    # 断言检查
    assert len(processed_net.edges()) == 3, "应生成3条新边"
    assert len(processed_net.nodes()) == 3, "应新增2个节点"
    assert any(data["length"] == 1.0 for u, v, data in processed_net.edges(data=True)), "每条新边长度应为3.0/3=1.0"
    assert not any(u == v for u, v in processed_net.edges()), "处理后不应存在自环边"


def test_get_matched_adjacencies():
    pass

def function_test():
    #处理重边的测试
    test_insert_one_node_into_duplicate_edges()
    #处理单边的测试
    test_change_single_edge_into_double()
    #处理自环的测试
    test_insert_two_nodes_into_selfloop()

    
    
#程序比较复杂，所以先确保函数功能没问题，再去验证逻辑正确
def test_metric_flip():
    #函数功能测试
    function_test()
    # 线性拓扑测试
    linear1 = pd.DataFrame(
        data=[["a", "b", 1.2, False],
              ["b", "c", 1.6, False]],
        columns=["from", "to", "length", "directed"]
    )
    linear2 = pd.DataFrame(
        data=[
            ["a", "b", 0.1, False],
            ["b", "c", 1.4, False],
            ["c", "d", 1.8, False]
        ],
        columns=["from", "to", "length", "directed"],
    )
    #分支拓扑逻辑逻辑验证
    bifurcating1=pd.DataFrame(data=[
        ["a",   "b", 1, True],
        ["b",   "c", 2, True],
        ["b",   "d", 3, True],
    ],
        columns=["from", "to", "length", "directed"])
    bifurcating2=pd.DataFrame(data=[
        ["b",   "a", 1, False],
        ["b",   "c", 2, False],
        ["b",   "d", 3, False],
        ["a",   "x", 4, False],
        ["c",   "y", 5, False],
        ["d",   "z", 6, False]

    ],
        columns=["from", "to", "length", "directed"])
    # 转换为 networkx 图
    def df_to_nx(df):
        G = nx.Graph()
        for _, row in df.iterrows():
            G.add_edge(row["from"], row["to"], length=row["length"])
            # 确保无向图的双向连接
            if not row["directed"]:
                G.add_edge(row["to"], row["from"], length=row["length"])
        return G
    # networkx接口
    # gr = nx.from_pandas_edgelist(net, source="from", target="to", create_using=nx.DiGraph if directed else nx.Graph)  # 构造graph
    # 转化为networkx
    linear1_to_networkx1 = df_to_nx(linear1)
    linear2_to_networkx2 = df_to_nx(linear2)
    bifurcating1_to_networkx1=df_to_nx(bifurcating1)
    bifurcating2_to_networkx2=df_to_nx(bifurcating2)
    # 调用函数测试返回得分
    #线性
    result1 = calculate_edge_flip(linear1_to_networkx1, linear2_to_networkx2, simplify=False)
    result2 = calculate_edge_flip(linear2_to_networkx2, linear1_to_networkx1, simplify=False)
    result3=calculate_edge_flip(linear1_to_networkx1, linear1_to_networkx1, simplify=False)
    #分支
    result4=calculate_edge_flip(bifurcating1_to_networkx1, bifurcating2_to_networkx2, simplify=False)
    result5=calculate_edge_flip( bifurcating2_to_networkx2, bifurcating1_to_networkx1, simplify=False)
    result6=calculate_edge_flip(bifurcating1_to_networkx1, bifurcating1_to_networkx1, simplify=False)
    #线性
    assert np.isclose(result1,4/5)
    assert np.isclose(result2,4/5)
    assert np.isclose(result3,1)
    #分支
    #assert np.isclose(result6,1)
    #循环

if __name__ == "__main__":
    pytest.main(["-v", __file__])
