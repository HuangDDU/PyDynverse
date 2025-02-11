import pytest
import pydynverse as pdv
import numpy as np
import pandas as pd

from ..wrap.test_wrap_add_dimred_projection import get_test_wrap_data


def test_project_to_segments():
    # 复用test_wrap_add_dimred_projection中的数据
    test_wrap_data = get_test_wrap_data()
    milestone_network = test_wrap_data["milestone_network"]
    dimred_milestones = test_wrap_data["dimred_milestones"]
    dimred = test_wrap_data["dimred"]
    expected_progressions = test_wrap_data["expected_progressions"]

    segment_start = dimred_milestones.loc[milestone_network["from"]]
    segment_end = dimred_milestones.loc[milestone_network["to"]]
    x = dimred

    # 预期输出
    expected_x_proj = pd.DataFrame(
        columns=["cell_id", "comp_1", "comp_2"],
        data=[
            ["a", 0.0, 1.0],
            ["b", 0.8, 1.0],
            ["c", 1.2, 1.0],
            ["d", 2.0, 1.0],
            ["e", 1.0, 1.5],
            ["f", 2.4, 1.0],
        ]
    ).set_index("cell_id", drop=True).values # 投影后的坐标
    expected_distance = np.array([0.25, 0.25, 0.25, 0.25, 0.04, 0.25]) # 源点到投影点的距离平方
    expected_segment = np.array([1, 1, 3, 3, 2, 4]) # 投影到在第几条边上
    expected_progression = np.array(expected_progressions["percentage"]) # 投影点在所在边额比例

    # 执行
    out = pdv.util.project_to_segments(x, segment_start, segment_end)

    assert (out["x_proj"] == expected_x_proj).any()
    assert (out["distance"] == expected_distance).any()
    assert (out["segment"] == expected_segment).any()
    assert (out["progression"] == expected_progression).any()


if __name__ == "__main__":
    pytest.main(["-v", __file__])
