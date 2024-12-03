import pytest
import pandas as pd
import pydynverse as pdv


def test_convert_milestone_percentages_to_progressions():
    id = "add_trajectory_new"
    cell_ids = ["a", "b", "c", "d", "e", "f"]
    milestone_ids = ["A", "B", "C", "D", "E", "F", "G"]
    milestone_network = pd.DataFrame(
        columns=["from", "to", "length", "directed"],
        data=[
            ["A", "B", 1, True],
            ["B", "C", 2, True],
            ["B", "D", 3, True],
            ["C", "E", 4, True],
            ["D", "F", 5, True],
            ["E", "G", 6, True],
            ["F", "G", 7, True]
        ]
    )
    milestone_percentages = pd.DataFrame(
        columns=["cell_id", "milestone_id", "percentage"],
        data=[
            ["a", "A", .8],
            ["a", "B", .2],
            ["b", "B", .3],
            ["b", "C", .2],
            ["b", "D", .5],
            ["c", "C", 0],
            ["c", "E", 1],
            ["d", "E", .5],
            ["d", "G", .5],
            ["e", "E", .9],
            ["e", "G", .1],
            ["f", "F", .6],
            ["f", "G", .4]
        ]
    )
    expected_progressions_percentage = [0.2, 0.2, 0.5, 1.0, 0.5, 0.1, 0.4]
    progressions = pdv.wrap.convert_milestone_percentages_to_progressions(
        cell_ids,
        milestone_ids,
        milestone_network,
        milestone_percentages
    )

    assert (progressions["percentage"]==expected_progressions_percentage).all(), "progressions percentage calculation failed"
    

if __name__ == "__main__":
    pytest.main()
    
