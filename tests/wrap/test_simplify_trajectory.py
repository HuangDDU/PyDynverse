import pytest
import pydynverse as pdv

import pandas as pd
from ..test_util import compare_dataframes_closely


def test_simplify_trajectory():
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

    pdv.wrap.simplify_trajectory(trajectory)

    # 预期输出
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

    assert trajectory["milestone_network"].equals(expected_milestone_network)
    assert compare_dataframes_closely(trajectory["progressions"], expected_progressions, on_columns="percentage")


if __name__ == "__main__":
    pytest.main(["-v", __file__])
