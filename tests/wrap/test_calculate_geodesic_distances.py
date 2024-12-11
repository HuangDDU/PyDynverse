import pytest

import pandas as pd
import pydynverse as pdv


def test_calculate_geodesic_distances():
    
    # 从test_wrap_add_waypoints.py导入测试样例数据
    from .test_wrap_add_waypoints import get_test_wrap_data
    dataset, milestone_network, divergence_regions, milestone_percentages = get_test_wrap_data()


if __name__ == "__main__":
    pytest.main(["-v", __file__])
