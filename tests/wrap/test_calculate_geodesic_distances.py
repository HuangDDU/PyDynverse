import pytest
import pydynverse as pdv

import time
from pydynverse.wrap.calculate_geodesic_distances import calculate_geodesic_distances

# def test_calculate_geodesic_distances():

#     # 从test_wrap_add_waypoints.py导入测试样例数据
#     from .test_wrap_add_waypoints import get_test_wrap_data
#     dataset, milestone_network, divergence_regions, milestone_percentages = get_test_wrap_data()


def test_calculate_geodesic_distances():
    dataset = pdv.data.load_simulation_data()
    out = calculate_geodesic_distances(dataset, graph_package="igraph")


if __name__ == "__main__":
    pytest.main(["-v", __file__])
