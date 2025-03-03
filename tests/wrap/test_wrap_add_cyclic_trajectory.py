import pandas as pd
import pytest
import pydynverse as pdv

from .test_wrap_add_linear_trajectory import get_test_wrap_data
from ..test_util import compare_dataframes_closely


def test_wrap_add_cyclic_trajectory():
    # input data
    test_wrap_data = get_test_wrap_data()
    dataset = test_wrap_data["dataset"]
    cell_ids = test_wrap_data["cell_ids"]
    pseudotime = test_wrap_data["pseudotime"]

    # execute function
    trajectory = pdv.wrap.add_cyclic_trajectory(
        dataset=dataset,
        pseudotime=pseudotime,
    )

    # expected result
    expected_milestone_ids = ["A", "B", "C"]
    expected_milestone_network = pd.DataFrame(
        columns=["from", "to", "length", "directed",],
        data=[
            ["A", "B", 1, False],
            ["B", "C", 1, False],
            ["C", "A", 1, False]
        ],
    )
    expected_progressions = pd.DataFrame(
        columns=["cell_id", "from", "to", "percentage"],
        data=[
            ["a", "A", "B", 0],
            ["b", "A", "B", 0.3],
            ["c", "B", "C", 0.2],
            ["d", "B", "C", 0.5],
            ["e", "C", "A", 0.4],
            ["f", 'C', "A", 1],
        ]
    )

    # assert
    assert trajectory["milestone_ids"] == expected_milestone_ids
    assert compare_dataframes_closely(trajectory["milestone_network"], expected_milestone_network, on_columns=["from", "to"])
    assert compare_dataframes_closely(trajectory["progressions"], expected_progressions, on_columns=["cell_id"])


if __name__ == "__main__":
    pytest.main(["-v", __file__])
