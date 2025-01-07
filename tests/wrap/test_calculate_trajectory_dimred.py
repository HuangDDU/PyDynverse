import pytest
import pydynverse as pdv

from .test_wrap_add_waypoints import get_test_wrap_data


def test_calculate_trajectory_dimred():
    # 复用test_wrap_add_waypoints数据
    dataset, milestone_network, divergence_regions, milestone_percentages = get_test_wrap_data()

    # 添加轨迹
    trajectory = pdv.wrap.add_trajectory(
        dataset,
        milestone_network=milestone_network,
        divergence_regions=divergence_regions,
        milestone_percentages=milestone_percentages
    )
    dimred_traj = pdv.wrap.calculate_trajectory_dimred(trajectory)
    assert dimred_traj["milestone_positions"].shape == (len(trajectory["milestone_ids"]), 2)
    assert dimred_traj["cell_positions"].shape == (len(trajectory["cell_ids"]), 2)
    assert True

if __name__ == "__main__":
    pytest.main(["-v", __file__])