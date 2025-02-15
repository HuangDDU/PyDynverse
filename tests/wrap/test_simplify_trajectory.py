import pytest
import pydynverse as pdv

import pandas as pd
from ..test_util import compare_dataframes_closely


def get_test_data():
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


def test_simplify_trajectory():

    test_data = get_test_data()
    trajectory = test_data["trajectory"]

    pdv.wrap.simplify_trajectory(trajectory)

    # 预期输出
    expected_milestone_network = test_data["expected_milestone_network"]
    expected_progressions = test_data["expected_progressions"]

    assert trajectory["milestone_network"].equals(expected_milestone_network)
    assert compare_dataframes_closely(trajectory["progressions"], expected_progressions, on_columns="percentage")


def test_simplify_trajectory_with_undirection():
    test_data = get_test_data()
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
    assert compare_dataframes_closely(trajectory["progressions"], expected_progressions, on_columns="percentage")


if __name__ == "__main__":
    pytest.main(["-v", __file__])
