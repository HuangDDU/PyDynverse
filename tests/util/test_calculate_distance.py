import pytest
from pydynverse.util import calculate_distance, list_similarity_methods
import pandas as pd


def test_calculate_distance():
    method_list = list_similarity_methods()
    # 对应../wrap/test_wrap_add_waypoints.py中的XYZ区域
    X = pd.DataFrame(
        data=[
            [0, 0, 0.2],
            [0, 0, 1],
            [0, 0.2, 0.5],
            [0, 0, 0],
            [0, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [0, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ],
        index=["c", "d", "e", "MILESTONE_BEGIN_WX_Y", "MILESTONE_BEGIN_WX_Z", "MILESTONE_END_WX_Y", "MILESTONE_END_WX_Z", "MILESTONE_X", "MILESTONE_Y", "MILESTONE_Z"],
        columns=["MILESTONE_X", "MILESTONE_Y", "MILESTONE_Z"]
    )

    wp_cells = ["MILESTONE_BEGIN_WX_Y", "MILESTONE_BEGIN_WX_Z", "MILESTONE_END_WX_Y", "MILESTONE_END_WX_Z"]
    tent = ["MILESTONE_X", "MILESTONE_Y", "MILESTONE_Z"]
    Y = X.loc[wp_cells+tent]

    # 计算距离
    distances = calculate_distance(X, Y, method="manhattan")
    
    # 预期构造结果
    expected_distances = pd.DataFrame(
        data=[
            [0.2, 0.2, 1.2, 0.8, 0.2, 1.2, 0.8],
            [1.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0],
            [0.7, 0.7, 1.3, 0.7, 0.7, 1.3, 0.7],
            [0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0],
            [0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0],
            [1.0, 1.0, 0.0, 2.0, 1.0, 0.0, 2.0],
            [1.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0],
            [0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0],
            [1.0, 1.0, 0.0, 2.0, 1.0, 0.0, 2.0],
            [1.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0],
        ],
        index=X.index,
        columns=Y.index
    )

    assert distances.equals(expected_distances)


if __name__ == "__main__":
    pytest.main(["-v", __file__])
