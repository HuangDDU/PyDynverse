import pytest
import pydynverse as pdv

import pandas as pd
import networkx as nx


def test_simplify_networkx_network():
    test = {
        "name": "directed_linear",
        "directed": True,
        "net": pd.DataFrame(
            data=[
                ["1", "2"],
                ["2", "3"],
                ["3", "4"],
            ],
            columns=["from", "to"]
        ),
        "expected_net": pd.DataFrame(
            data=[["1", "4", 3, True]],
            columns=["from", "to", "weight", "directed"]
        )
    }
    net = test["net"]
    gr = nx.from_pandas_edgelist(net, source="from", target="to", create_using=nx.DiGraph() if test["directed"] else nx.Graph)  # 构造graph
    expected_net = test["expected_net"]  # 预期输出

    new_gr = pdv.wrap.simplify_networkx_network(gr)

    # 提取输出
    new_net = pd.DataFrame(new_gr.edges(data=True), columns=["from", "to", "attributes"])
    new_net = pd.concat([new_net.drop(columns=['attributes']), new_net["attributes"].apply(pd.Series)], axis=1)

    assert new_net.equals(expected_net)


if __name__ == "__main__":
    pytest.main(["-v", __file__])
