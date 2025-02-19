import pytest
import numpy as np
import pydynverse as pdv
import pandas as pd
from sklearn.metrics import jaccard_score
from pandas import DataFrame, Series

from ..wrap.test_wrap_add_waypoints import get_test_wrap_data


def create_test_data():
    dataset, milestone_network, divergence_regions, milestone_percentages = get_test_wrap_data()

    # 添加轨迹
    trajectory_ref = pdv.wrap.add_trajectory(
        dataset,
        milestone_network=milestone_network,
        divergence_regions=divergence_regions,
        milestone_percentages=milestone_percentages
    )

    milestone_percentages

    trajectory_predict = None

def test_metric_mapping():
    # 获取测试数据
    dataset, prediction = create_test_data()
    # 测试按里程碑映射
    result_milestones = pdv.eval.calculate_mapping_milestones(dataset, prediction, simplify=True)
    # 测试按照分支映射
    result_branches = pdv.eval.calculate_mapping_branches(dataset, prediction, simplify=True)



if __name__ == "__main__":
    pytest.main(["-v", __file__])
