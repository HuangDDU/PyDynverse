import pytest
from ..wrap.test_wrap_add_waypoints import get_test_wrap_data

def test_plot_dimred():
    # 测试样例数据
    dataset, milestone_network, divergence_regions, milestone_percentages = get_test_wrap_data()


if __name__ == "__main__":
    pytest.main(["-v", __file__])
