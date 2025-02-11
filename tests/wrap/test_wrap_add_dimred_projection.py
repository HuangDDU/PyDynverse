import pytest
import pydynverse as pdv
import pandas as pd

from .test_wrap_add_waypoints import get_test_wrap_data as get_test_wrap_data_ref


def get_test_wrap_data():
    # 复用test_wrap_add_waypoints中的数据
    dataset, milestone_network, divergence_regions, milestone_percentages = get_test_wrap_data_ref()
    # expected_milestone_percentages = milestone_percentages[milestone_percentages["cell_id"] == "e"]
    expected_milestone_percentages = pd.DataFrame(
        columns=["cell_id", "milestone_id", "percentage"],
        data=[
            ["e", "X", 0.5],
            ["e", "Y", 0.5]
        ]
    )
    expected_milestone_percentages = pd.concat([expected_milestone_percentages, milestone_percentages.query("cell_id != 'e'")])
    expected_progressions = pd.DataFrame(
        columns=["cell_id", "from", "to", "percentage"],
        data=[
            ["a", "W", "X", 0],
            ["b", "W", "X", 0.8],
            ["c", "X", "Z", 0.2],
            ["d", "X", "Z", 1],
            ["e", "X", "Y", 0.5],
            ["f", "Z", "A", 0.2],
        ]
    )

    # 手动指定2维降维
    dimred = pd.DataFrame(
        columns=["cell_id", "comp_1", "comp_2"],
        data=[
            ["a", 0, 1.5],
            ["b", 0.8, 0.5],
            ["c", 1.2, 0.5],
            ["d", 2, 0.5],
            ["e", 1.2, 1.5],
            ["f", 2.4, 1.5],
        ]
    )
    dimred_milestones = pd.DataFrame(
        columns=["milestone_id", "comp_1", "comp_2"],
        data=[
            ["W", 0, 1],
            ["X", 1, 1],
            ["Y", 1, 2],
            ["Z", 2, 1],
            ["A", 4, 1],
        ]
    )
    dimred.set_index("cell_id", inplace=True)
    dimred_milestones.set_index("milestone_id", inplace=True)

    test_wrap_data = {
        "dataset": dataset,
        "milestone_network": milestone_network,
        "divergence_regions": divergence_regions,
        "milestone_percentages": milestone_percentages,
        "dimred": dimred,
        "dimred_milestones": dimred_milestones,
        "expected_milestone_percentages": expected_milestone_percentages,
        "expected_progressions": expected_progressions,
    }

    return test_wrap_data


def test_add_dimred_projection():
    test_wrap_data = get_test_wrap_data()
    dataset = test_wrap_data["dataset"]
    milestone_network = test_wrap_data["milestone_network"]
    dimred = test_wrap_data["dimred"]
    dimred_milestones = test_wrap_data["dimred_milestones"]
    expected_milestone_percentages = test_wrap_data["expected_milestone_percentages"]
    expected_progression = test_wrap_data["expected_progression"]

    trajectory = pdv.wrap.add_dimred_projection(
        dataset=dataset,
        milestone_network=milestone_network,
        dimred=dimred,
        dimred_milestones=dimred_milestones
    )
    # assert trajectory["milestone_percentages"].equals(expected_milestone_percentages)
    # assert trajectory["progression"].equals(expected_progression)


if __name__ == "__main__":
    pytest.main(["-v", __file__])
