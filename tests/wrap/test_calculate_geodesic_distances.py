import pytest
import pydynverse as pdv

import time
import pandas as pd
from pydynverse.wrap.calculate_geodesic_distances import calculate_geodesic_distances

# def test_calculate_geodesic_distances():

#     # 从test_wrap_add_waypoints.py导入测试样例数据
#     from .test_wrap_add_waypoints import get_test_wrap_data
#     dataset, milestone_network, divergence_regions, milestone_percentages = get_test_wrap_data()

def get_data():
    id = "test_calculate_geodesic_distances"
    cell_ids = ["a", "b", "c", "d", "e"]
    milestone_ids = ["W", "X", "Y", "Z"]
    milestone_network = pd.DataFrame(
        columns=["from", "to", "length", "directed"],
        data=[
            ["W", "X", 2, True],
            ["X", "Y", 3, True],
            ["X", "Z", 4, True],
        ]
    )
    divergence_regions = pd.DataFrame(
        columns=["divergence_id", "milestone_id", "is_start"],
        data=[
            ["XYZ", "X", True],
            ["XYZ", "Y", False],
            ["XYZ", "Z", False]
        ]
    )
    milestone_percentages = pd.DataFrame(
        columns=["cell_id", "milestone_id", "percentage"],
        data=[
            ["a", "W", 0.9],
            ["a","X",0.1],
            ["b", "W", 0.2],
            ["b", "X", 0.8],
            ["c", "X", 0.8],
            ["c", "Z", 0.2],
            ["d", "Z", 0.1],
            ["d", "X", 0.2],
            ["d", "Y", 0.7],
            ["e", "X", 0.3],
            ["e", "Y", 0.2],
            ["e", "Z", 0.5],
        ]
    )
    dataset = pdv.wrap.wrap_data(id=id, cell_ids=cell_ids)
    test_data = {
        "dataset": dataset,
        "milestone_network": milestone_network,
        "divergence_regions": divergence_regions,
        "milestone_percentages": milestone_percentages,
    }

    dataset=test_data["dataset"]
    milestone_network=test_data["milestone_network"]
    divergence_regions=test_data["divergence_regions"]
    milestone_percentages=test_data["milestone_percentages"]
    # 添加轨迹
    test_trajectory = pdv.wrap.add_trajectory(
        dataset,
        milestone_network=milestone_network,
        divergence_regions=divergence_regions,
        milestone_percentages=milestone_percentages
    )
    return test_trajectory

def test_calculate_geodesic_distances():
    # dataset = pdv.data.load_simulation_data()
    trajectory = get_data()
    out = calculate_geodesic_distances(trajectory, graph_package="igraph")

    assert True 


if __name__ == "__main__":
    pytest.main(["-v", __file__])
