import pytest
import pydynverse as pdv

import pandas as pd
import networkx as nx
from ..test_util import compare_dataframes_closely


def test_simplify_networkx_network():
    directed = True
    net = pd.DataFrame(
        data=[
            ["A", "B"],
            ["B", "C"],
            ["C", "D"],
        ],
        columns=["from", "to"]
    )
    edge_points = pd.DataFrame(
        data=[
            ["a", "A", "B", 0.3],
            ["b", "A", "B", 0.6],
            ["c", "B", "C", 0.2],
            ["d", "B", "C", 0.8],
            ["e", "C", "D", 0.4],
        ],
        columns=["id", "from", "to", "percentage"]
    )
    gr = nx.from_pandas_edgelist(net, source="from", target="to", create_using=nx.DiGraph() if directed else nx.Graph)  # 构造graph

    out = pdv.wrap.simplify_networkx_network(gr, edge_points=edge_points)
    new_gr = out["gr"]
    new_edge_points = out["edge_points"]

    # 提取输出
    new_net = pd.DataFrame(new_gr.edges(data=True), columns=["from", "to", "attributes"])
    new_net = pd.concat([new_net.drop(columns=['attributes']), new_net["attributes"].apply(pd.Series)], axis=1)

    # 预期输出    
    expected_net = pd.DataFrame(
        data=[["A", "D", 3, True]],
        columns=["from", "to", "weight", "directed"]
    )
    expected_edge_points = pd.DataFrame(
        data=[
            ["a", "A", "D", 0.1],
            ["b", "A", "D", 0.2],
            ["c", "A", "D", 0.4],
            ["d", "A", "D", 0.6],
            ["e", "A", "D", 0.8],
        ],
        columns=["id", "from", "to", "percentage"]
    )

    assert new_net.equals(expected_net)
    assert compare_dataframes_closely(new_edge_points, expected_edge_points, on_columns="percentage")



if __name__ == "__main__":
    pytest.main(["-v", __file__])
