import pytest
import pandas as pd
import pydynverse as pdv


# 数据部分单独提取出来，其他测试用例还能用，例如test_convert_milestone_percentages_to_progressions.py
def get_test_wrap_data():
    id = "add_trajectory_new"
    cell_ids = ["a", "b", "c", "d", "e", "f"]
    milestone_ids = ["A", "B", "C", "D", "E", "F", "G"]
    milestone_network = pd.DataFrame(
        columns=["from", "to", "length", "directed"],
        data=[
            ["A", "B", 1, True],
            ["B", "C", 2, True],
            ["B", "D", 3, True],
            ["C", "E", 4, True],
            ["D", "F", 5, True],
            ["E", "G", 6, True],
            ["F", "G", 7, True]
        ]
    )
    milestone_percentages = pd.DataFrame(
        columns=["cell_id", "milestone_id", "percentage"],
        data=[
            ["a", "A", .8],
            ["a", "B", .2],
            ["b", "B", .3],
            ["b", "C", .2],
            ["b", "D", .5],
            ["c", "C", 0],
            ["c", "E", 1],
            ["d", "E", .5],
            ["d", "G", .5],
            ["e", "E", .9],
            ["e", "G", .1],
            ["f", "F", .6],
            ["f", "G", .4]
        ]
    )
    divergence_regions = pd.DataFrame(
        columns=["divergence_id", "milestone_id", "is_start"],
        data=[
            ["BCD", "B", True],
            ["BCD", "C", False],
            ["BCD", "D", False]
        ]
    )

    dataset = pdv.wrap.wrap_data(id=id, cell_ids=cell_ids)
    return dataset, milestone_network, milestone_percentages, divergence_regions


def test_add_trajectory():
    # 测试样例数据
    dataset, milestone_network, milestone_percentages, divergence_regions = get_test_wrap_data()

    # 添加轨迹
    trajectory = pdv.wrap.add_trajectory(
        dataset,
        milestone_network=milestone_network,
        divergence_regions=divergence_regions,
        milestone_percentages=milestone_percentages
    )

    # 这里同时测试了gather_cells_at_milestones函数
    gathered_trajectory = pdv.wrap.gather_cells_at_milestones(trajectory)
    # 聚合后的milestone_percentages DataFrame测试
    gathered_milestone_percentages = gathered_trajectory["milestone_percentages"]
    expected_milestone_id_list = ["A", "D", "E", "E", "E", "F"]
    assert (gathered_milestone_percentages["percentage"] == 1).all(), "gathered_milestone_percentages percentage error"
    assert (gathered_milestone_percentages["milestone_id"] == expected_milestone_id_list).all(), "gathered_milestone_percentages milestone_id error "
    # assert type(milestone_network) == pd.DataFrame, "类型错误"


if __name__ == "__main__":
    pytest.main(["-v", __file__])
