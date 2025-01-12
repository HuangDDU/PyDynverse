import pytest
import pydynverse as pdv
from .test_calculate_trajectory_dimred import get_test_wrap_data


def test_group_onto_trajectory_edges():
    trajectory = get_test_wrap_data()
    group = pdv.wrap.group_onto_trajectory_edges(trajectory)
    # excepted_group = ["W->X", "W->X", "X->Z", "X->Z", "X->Z", "Z->A"]
    # TODO: wrap.convert_milestone_percentages_to_progressions对于milestone上的细胞的progression: X->X
    excepted_group = ["W->W", "W->X", "X->Z", "Z->Z", "X->Z", "Z->A"]
    assert group == excepted_group


def test_group_onto_nearest_milestones():
    trajectory = get_test_wrap_data()
    group = pdv.wrap.group_onto_nearest_milestones(trajectory)
    excepted_group = ["W", "X", "X", "Z", "Z", "Z"]
    assert group == excepted_group


if __name__ == "__main__":
    pytest.main(["-v", __file__])
