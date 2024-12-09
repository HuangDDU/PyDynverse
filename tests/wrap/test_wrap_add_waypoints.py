import pytest

import numpy as np
import pandas as pd
import pydynverse as pdv


def get_test_wrap_data():
    id = "wrap_add_waypoints"
    cell_ids = ["a", "b", "c", "d", "e", "f"]
    milestone_ids = ["W", "X", "Y", "Z", "A"]

    milestone_network = pd.DataFrame(
        columns=["from", "to", "length", "directed"],
        data=[
            ["W", "X", 1, True],
            ["X", "Y", 1, True],
            ["X", "Z", 1, True],
            ["Z", "A", 2, True]
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
            ["a", "W", 1],
            ["b", "W", .2],
            ["b", "X", .8],
            ["c", "X", .8],
            ["c", "Z", .2],
            ["d", "Z", 1],
            ["e", "X", .3],
            ["e", "Y", .2],
            ["e", "Z", .5],
            ["f", "Z", .8],
            ["f", "A", .2],
        ]
    )

    dataset = pdv.wrap.wrap_data(id=id, cell_ids=cell_ids)
    return dataset, milestone_network, divergence_regions, milestone_percentages


def test_wrap_add_waypoints():
    # 测试样例数据
    dataset, milestone_network, divergence_regions, milestone_percentages = get_test_wrap_data()

    # 添加轨迹
    trajectory = pdv.wrap.add_trajectory(
        dataset,
        milestone_network=milestone_network,
        divergence_regions=divergence_regions,
        milestone_percentages=milestone_percentages
    )

    # 选择waypoints
    wp = pdv.wrap.select_waypoints(
        trajectory=trajectory,
        trafo=lambda x: x,
        resolution=1
    )  # milstone边上每间隔1个单位找一个waypoint

    # 预期构造结果
    expected_milestone_percentages = pd.DataFrame(
        columns=["waypoint_id", "milestone_id", "percentage"],
        data=[
            ["MILESTONE_BEGIN_WW_X", "W", 1],
            ["MILESTONE_BEGIN_WW_X", "X", 0],
            ["MILESTONE_BEGIN_WX_Y", "X", 1],
            ["MILESTONE_BEGIN_WX_Y", "Y", 0],
            ["MILESTONE_BEGIN_WX_Z", "X", 1],
            ["MILESTONE_BEGIN_WX_Z", "Z", 0],
            ["MILESTONE_BEGIN_WZ_A", "Z", 1],
            ["MILESTONE_BEGIN_WZ_A", "A", 0],
            ["MILESTONE_END_WW_X", "W", 0],
            ["MILESTONE_END_WW_X", "X", 1],
            ["MILESTONE_END_WX_Y", "X", 0],
            ["MILESTONE_END_WX_Y", "Y", 1],
            ["MILESTONE_END_WX_Z", "X", 0],
            ["MILESTONE_END_WX_Z", "Z", 1],
            ["MILESTONE_END_WZ_A", "Z", 0],
            ["MILESTONE_END_WZ_A", "A", 1],
            ["W8", "Z", 0.5],
            ["W8", "A", 0.5],
        ]
    )
    expected_progression = pd.DataFrame(
        columns=["from", "to", "percentage", "waypoint_id"],
        data=[
            ["W", "X", 0, "MILESTONE_BEGIN_WW_X"],
            ["W", "X", 1, "MILESTONE_END_WW_X"],
            ["X", "Y", 0, "MILESTONE_BEGIN_WX_Y"],
            ["X", "Y", 1, "MILESTONE_END_WX_Y"],
            ["X", "Z", 0, "MILESTONE_BEGIN_WX_Z"],
            ["X", "Z", 1, "MILESTONE_END_WX_Z"],
            ["Z", "A", 0, "MILESTONE_BEGIN_WZ_A"],
            ["Z", "A", 0.5, "W8"],
            ["Z", "A", 1, "MILESTONE_END_WZ_A"],
        ]
    )
    expected_geodesic_distances = pd.DataFrame(
        columns=dataset["cell_ids"],
        index=["MILESTONE_BEGIN_WW_X",
               "MILESTONE_END_WW_X",
               "MILESTONE_BEGIN_WX_Y",
               "MILESTONE_END_WX_Y",
               "MILESTONE_BEGIN_WX_Z",
               "MILESTONE_END_WX_Z",
               "MILESTONE_BEGIN_WZ_A",
               "W8",
               "MILESTONE_END_WZ_A",],
        data=[
            [0.0, 0.8, 1.2, 2.0, 1.7, 2.4],
            [1.0, 0.2, 0.2, 1.0, 0.7, 1.4],
            [1.0, 0.2, 0.2, 1.0, 0.7, 1.4],
            [2.0, 1.2, 1.2, 2.0, 1.3, 2.4],
            [1.0, 0.2, 0.2, 1.0, 0.7, 1.4],
            [2.0, 1.2, 0.8, 0.0, 0.7, 0.4],
            [2.0, 1.2, 0.8, 0.0, 0.7, 0.4],
            [3.0, 2.2, 1.8, 1.0, 1.7, 0.6],
            [4.0, 3.2, 2.8, 2.0, 2.7, 1.6],
        ]
    )
    expected_waypoint_network = pd.DataFrame(
        columns=["from", "to", "from_milestone_id", "to_milestone_id"],
        data=[
            ["MILESTONE_BEGIN_WW_X", "MILESTONE_END_WW_X", "W", "X"],
            ["MILESTONE_BEGIN_WX_Y", "MILESTONE_END_WX_Y", "X", "Y"],
            ["MILESTONE_BEGIN_WX_Z", "MILESTONE_END_WX_Z", "X", "Z"],
            ["MILESTONE_BEGIN_WZ_A", "W8", "Z", "A"],
            ["W8", "MILESTONE_END_WZ_A", "Z", "A"],
        ]
    )
    expected_waypoints = pd.DataFrame(
        columns=["waypoint_id", "milestone_id"],
        data=[
            ["MILESTONE_BEGIN_WW_X", "W"],
            ["MILESTONE_BEGIN_WX_Y", "X"],
            ["MILESTONE_BEGIN_WX_Z", "X"],
            ["MILESTONE_BEGIN_WZ_A", "Z"],
            ["MILESTONE_END_WW_X", "X"],
            ["MILESTONE_END_WX_Y", "Y"],
            ["MILESTONE_END_WX_Z", "Z"],
            ["MILESTONE_END_WZ_A", "A"],
            ["W8", None],
        ]
    )

    # 开始断言
    # 关于milestone_percentages, 两个dataframe行顺序不一致, 需要通过waypoint_id和milestone_id联合排序后再比较
    sorted_expected_milestone_percentages = expected_milestone_percentages\
        .sort_values(by=["waypoint_id", "milestone_id"])\
        .reset_index(drop=True)
    sorted_wp_milestone_percentages = wp["milestone_percentages"]\
        .rename(columns={"cell_id": "waypoint_id"})\
        .sort_values(by=["waypoint_id", "milestone_id"])\
        .reset_index(drop=True)
    assert sorted_wp_milestone_percentages.equals(sorted_expected_milestone_percentages)

    assert wp["progressions"].equals(expected_progression)

    # assert wp["geodesic_distances"].equals(expected_geodesic_distances) # 这里调试看的数值是完全一样的，函数调用计算距离经过多次运算后的精度不一致
    # expected_geodesic_distances.values[1,1]=0.2, wp["geodesic_distances"].values[1,1]=0.19999999999999996
    assert np.allclose(wp["geodesic_distances"].values, expected_geodesic_distances.values) # 误差在可容忍的范围内

    # waypoint_network关于from与to联合排序后再比较
    sorted_expected_waypoint_network = expected_waypoint_network\
        .sort_values(by=["from", "to"])\
        .reset_index(drop=True)
    sorted_wp_waypoint_network = wp["waypoint_network"]\
        .sort_values(by=["from", "to"])\
        .reset_index(drop=True)
    assert sorted_wp_waypoint_network.equals(sorted_expected_waypoint_network)

    # waypoints关waypoint_id排序后再比较
    sorted_expected_waypoints = expected_waypoints\
        .sort_values(by=["waypoint_id"])\
        .reset_index(drop=True)
    sorted_wp_waypoints = wp["waypoints"]\
        .sort_values(by=["waypoint_id"])\
        .reset_index(drop=True)
    assert sorted_wp_waypoints.equals(sorted_expected_waypoints)


if __name__ == "__main__":
    pytest.main(["-v", __file__])
