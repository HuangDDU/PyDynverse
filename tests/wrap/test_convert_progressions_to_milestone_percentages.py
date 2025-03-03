import pytest
import pandas as pd
import pydynverse as pdv

from ..test_util import compare_dataframes_closely


def test_convert_progressions_to_milestone_percentages():
    cell_ids = ["a", "b", "c"]
    milestone_ids = ["milestone_begin", "A", "B", "C"]
    milestone_network = pd.DataFrame(
        columns=["from", "to", "length", "directed"],
        data=[
            ["milestone_begin", "A", 1, True],
            ["milestone_begin", "B", 1, True],
            ["milestone_begin", "C", 1, True],
        ],
    )
    end_state_probabilities = pd.DataFrame(
        columns=["cell_id", "A", "B", "C"],
        data=[
            ["a", 0.5, 0.2, 0.2],
            ["b", 0.2, 0.5, 0.2],
            ["c", 0.2, 0.2, 0.5],
        ],
    )
    progressions = end_state_probabilities.melt(id_vars=["cell_id"], var_name="to", value_name="percentage")
    progressions["from"] = "milestone_begin"

    milestone_percentages = pdv.wrap.convert_progressions_to_milestone_percentages(
        cell_ids=cell_ids,
        milestone_ids=milestone_ids,
        milestone_network=milestone_network,
        progressions=progressions,
    )

    expected_milestone_percentages = pd.DataFrame(
        columns=["cell_id", "milestone_id", "percentage"],
        data=[
            ["a", "milestone_begin", 0.1],  # for start milestone， percentage = 1 - sum(other end milestone percentages)
            ["a", "A", 0.5],
            ["a", "B", 0.2],
            ["a", "C", 0.2],
            ["b", "milestone_begin", 0.1],
            ["b", "A", 0.2],
            ["b", "B", 0.5],
            ["b", "C", 0.2],
            ["c", "milestone_begin", 0.1],
            ["c", "A", 0.2],
            ["c", "B", 0.2],
            ["c", "C", 0.5],
        ],
    )

    assert compare_dataframes_closely(milestone_percentages, expected_milestone_percentages, on_columns=["cell_id", "milestone_id"])


if __name__ == "__main__":
    pytest.main(["-v", __file__])
