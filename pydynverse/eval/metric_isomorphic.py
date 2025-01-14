import networkx as nx


def calc_isomorphic(net1, net2):
    # 图同构
    graph1 = nx.from_pandas_edgelist(net1, source="from", target="to")
    graph2 = nx.from_pandas_edgelist(net2, source="from", target="to")
    if nx.is_isomorphic(graph1, graph2):
        return 1
    else:
        return 0
