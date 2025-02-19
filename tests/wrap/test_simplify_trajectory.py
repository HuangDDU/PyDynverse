import pytest
import pydynverse as pdv

import pandas as pd
from ..test_util import compare_dataframes_closely


def get_test_data_linear():
    id = "directed_linear"
    cell_ids = ["a", "b", "c", "d", "e"]
    milestone_ids = ["A", "B", "C", "D"]
    milestone_network = pd.DataFrame(
        data=[
            ["A", "B", 1, True],
            ["B", "C", 1, True],
            ["C", "D", 1, True]
        ],
        columns=["from", "to", "length", "directed"],
    )
    progressions = pd.DataFrame(
        data=[
            ["a", "A", "B", 0.3],
            ["b", "A", "B", 0.6],
            ["c", "B", "C", 0.2],
            ["d", "B", "C", 0.8],
            ["e", "C", "D", 0.4],
        ],
        columns=["cell_id", "from", "to", "percentage"]
    )

    dataset = pdv.wrap.wrap_data(id=id, cell_ids=cell_ids)
    trajectory = pdv.wrap.add_trajectory(
        dataset=dataset,
        milestone_ids=milestone_ids,
        milestone_network=milestone_network,
        progressions=progressions
    )

    expected_milestone_network = pd.DataFrame(
        data=[["A", "D", 3, True]],
        columns=["from", "to", "length", "directed"],
    )
    expected_progressions = pd.DataFrame(
        data=[
            ["a", "A", "D", 0.1],
            ["b", "A", "D", 0.2],
            ["c", "A", "D", 0.4],
            ["d", "A", "D", 0.6],
            ["e", "A", "D", 0.8],
        ],
        columns=["cell_id", "from", "to", "percentage"]
    )

    test_data = {
        "id": id,
        "cell_ids": cell_ids,
        "milestone_ids": milestone_ids,
        "milestone_network": milestone_network,
        "progressions": progressions,
        "trajectory": trajectory,
        "expected_milestone_network": expected_milestone_network,
        "expected_progressions": expected_progressions,
    }
    return test_data


def test_simplify_trajectory_directed_linear():

    test_data = get_test_data_linear()
    trajectory = test_data["trajectory"]

    pdv.wrap.simplify_trajectory(trajectory)

    # 预期输出
    expected_milestone_network = test_data["expected_milestone_network"]
    expected_progressions = test_data["expected_progressions"]

    assert trajectory["milestone_network"].equals(expected_milestone_network)
    assert compare_dataframes_closely(trajectory["progressions"], expected_progressions, on_columns="cell_id")


def test_simplify_trajectory_undirected_linear():
    test_data = get_test_data_linear()
    id = test_data["id"]
    cell_ids = test_data["cell_ids"]
    milestone_ids = test_data["milestone_ids"]
    milestone_network = test_data["milestone_network"]
    progressions = test_data["progressions"]
    milestone_network["directed"] = False  # undirected graph
    dataset = pdv.wrap.wrap_data(id=id, cell_ids=cell_ids)
    trajectory = pdv.wrap.add_trajectory(
        dataset=dataset,
        milestone_ids=milestone_ids,
        milestone_network=milestone_network,
        progressions=progressions
    )

    pdv.wrap.simplify_trajectory(trajectory)

    # 预期输出
    expected_milestone_network = test_data["expected_milestone_network"]
    expected_milestone_network["directed"] = False
    expected_progressions = test_data["expected_progressions"]

    assert trajectory["milestone_network"].equals(expected_milestone_network)
    assert compare_dataframes_closely(trajectory["progressions"], expected_progressions, on_columns="cell_id")


def get_test_data_bifurcation():
    id = "directed_bifurcation"
    cell_ids = ["a", "b", "c", "d", "e", "f"]
    milestone_ids = ["A", "B", "C", "D", "E", "F", "G"]
    milestone_network = pd.DataFrame(
        data=[
            ["A", "B", 4, True],
            ["A", "C", 4, True],
            ["B", "D", 1, True],
            ["C", "E", 1, True],
            ["E", "F", 1, True],
            ["E", "G", 1, True],
        ],
        columns=["from", "to", "length", "directed"],
    )
    progressions = pd.DataFrame(
        data=[
            ["a", "A", "B", 0.5],
            ["b", "A", "C", 0.5],
            ["c", "B", "D", 0.5],
            ["d", "C", "E", 0.5],
            ["e", "E", "F", 0.5],
            ["f", "E", "G", 0.5],
        ],
        columns=["cell_id", "from", "to", "percentage"]
    )

    dataset = pdv.wrap.wrap_data(id=id, cell_ids=cell_ids)
    trajectory = pdv.wrap.add_trajectory(
        dataset=dataset,
        milestone_ids=milestone_ids,
        milestone_network=milestone_network,
        progressions=progressions
    )

    expected_milestone_network = pd.DataFrame(
        data=[
            ["A", "D", 5, True],
            ["A", "E", 5, True],
            ["E", "F", 1, True],
            ["E", "G", 1, True]
        ],
        columns=["from", "to", "length", "directed"],
    )
    expected_progressions = pd.DataFrame(
        data=[
            ["a", "A", "D", 0.4],
            ["b", "A", "E", 0.4],
            ["c", "A", "D", 0.9],
            ["d", "A", "E", 0.9],
            ["e", "E", "F", 0.5],
            ["f", "E", "G", 0.5],
        ],
        columns=["cell_id", "from", "to", "percentage"]
    )

    test_data = {
        "id": id,
        "cell_ids": cell_ids,
        "milestone_ids": milestone_ids,
        "milestone_network": milestone_network,
        "progressions": progressions,
        "trajectory": trajectory,
        "expected_milestone_network": expected_milestone_network,
        "expected_progressions": expected_progressions,
    }
    return test_data


def test_simplify_trajectory_directed_bifurcation():
    test_data = get_test_data_bifurcation()
    trajectory = test_data["trajectory"]

    pdv.wrap.simplify_trajectory(trajectory)

    # 预期输出
    expected_milestone_network = test_data["expected_milestone_network"]
    expected_progressions = test_data["expected_progressions"]

    assert trajectory["milestone_network"].equals(expected_milestone_network)
    assert compare_dataframes_closely(trajectory["progressions"], expected_progressions, on_columns="cell_id")


def test_simplify_trajectory_undirected_bifurcation():
    test_data = get_test_data_bifurcation()
    id = test_data["id"]
    cell_ids = test_data["cell_ids"]
    milestone_ids = test_data["milestone_ids"]
    milestone_network = test_data["milestone_network"]
    progressions = test_data["progressions"]
    milestone_network["directed"] = False  # undirected graph
    dataset = pdv.wrap.wrap_data(id=id, cell_ids=cell_ids)
    trajectory = pdv.wrap.add_trajectory(
        dataset=dataset,
        milestone_ids=milestone_ids,
        milestone_network=milestone_network,
        progressions=progressions
    )

    # 执行
    pdv.wrap.simplify_trajectory(trajectory)

    # 预期输出
    expected_milestone_network = pd.DataFrame(
        data=[
            ["D", "E", 10, False],
            ["E", "F", 1, False],
            ["E", "G", 1, False],
        ],
        columns=["from", "to", "length", "directed"],
    )
    expected_progressions = pd.DataFrame(
        data=[
            ["a", "D", "E", 0.3],
            ["b", "D", "E", 0.7],
            ["c", "D", "E", 0.05],
            ["d", "D", "E", 0.95],
            ["e", "E", "F", 0.5],
            ["f", "E", "G", 0.5],
        ],
        columns=["cell_id", "from", "to", "percentage"]
    )
    
    # assert
    assert trajectory["milestone_network"].equals(expected_milestone_network)
    assert compare_dataframes_closely(trajectory["progressions"], expected_progressions, on_columns="cell_id") # TODO: 这里暂时有问题，progression里出现了milestone_network中没有的milestone


if __name__ == "__main__":
    pytest.main(["-v", __file__])
