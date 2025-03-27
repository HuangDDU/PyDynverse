import pytest
import pydynverse as pdv

import numpy as np
import pandas as pd


def get_test_wrap_data():
    from tests.wrap.test_wrap_add_waypoints import get_test_wrap_data as get_test_wrap_data_ref
    # 这里复用了test_wrap_add_waypoint的数据
    test_wrap_data= get_test_wrap_data_ref()
    dataset = test_wrap_data["dataset"]
    milestone_network = test_wrap_data["milestone_network"]
    divergence_regions = test_wrap_data["divergence_regions"]
    milestone_percentages = test_wrap_data["milestone_percentages"]

    # 添加轨迹
    trajectory_ref = pdv.wrap.add_trajectory(
        dataset,
        milestone_network=milestone_network,
        divergence_regions=divergence_regions,
        milestone_percentages=milestone_percentages
    )

    new_mp1 = milestone_percentages.query("`cell_id`!='e'").copy()
    new_mp2 = pd.DataFrame(
        columns=["cell_id", "milestone_id", "percentage"],
        data=[
            ["e", "X", 0.2],
            ["e", "Y", 0.8],
        ]
    )
    new_milestone_percentages = pd.concat([new_mp1, new_mp2])
    trajectory_pre = pdv.wrap.add_trajectory(
        dataset.copy(),
        milestone_network=milestone_network,
        divergence_regions=divergence_regions,
        milestone_percentages=new_milestone_percentages
    )

    test_wrap_data = {
        "trajectory_ref": trajectory_ref,
        "trajectory_pre": trajectory_pre,
    }

    return test_wrap_data


def test_metric_mapping():

    # 1. 提取参考、预测轨迹
    test_wrap_data = get_test_wrap_data()
    trajectory_ref = test_wrap_data["trajectory_ref"]
    trajectory_pre = test_wrap_data["trajectory_pre"]

    # 2.调用接口获得模型结果
    result_milestones = pdv.eval.calculate_mapping_milestones(trajectory_ref, trajectory_pre)
    result_branches = pdv.eval.calculate_mapping_branches(trajectory_ref, trajectory_pre)

    # 3.期望的结果
    expected_milestones_recovery = 8/9
    expected_milestones_relevance = 3/4
    expected_milestones_F1 = 48/59
    # TODO: 手动计算后再验证
    
    expected_branches_F1 = None

    # 4.测试断言
    assert np.isclose(result_milestones["recovery_milestones"], expected_milestones_recovery)
    assert np.isclose(result_milestones["relevance_milestones"], expected_milestones_relevance)
    assert np.isclose(result_milestones["F1_milestones"], expected_milestones_F1)
    # assert   result_branches["F1_branches"]<1


if __name__ == "__main__":
    pytest.main(["-v", __file__])
