import pytest
import pandas as pd
import pydynverse as pdv


def test_convert_milestone_percentages_to_progressions():
    
    # 从test_wrap_add_trajectory.py提取测试样例数据
    from .test_wrap_add_trajectory import get_test_wrap_data
    dataset, milestone_network, milestone_percentages, divergence_regions = get_test_wrap_data()
      
    # 执行转换
    progressions = pdv.wrap.convert_milestone_percentages_to_progressions(
        cell_ids=dataset["cell_ids"],
        milestone_ids="",
        milestone_network=milestone_network,
        milestone_percentages=milestone_percentages
    )

    # 预期转换后的progressions_percentage
    expected_progressions_percentage = [0.2, 0.2, 0.5, 1.0, 0.5, 0.1, 0.4]

    # 判断转换后的结果与预期结果是否一致
    assert (progressions["percentage"]==expected_progressions_percentage).all(), "progressions percentage calculation failed"
    

if __name__ == "__main__":
    pytest.main(["-v", __file__])
