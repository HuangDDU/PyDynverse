import pytest
import pandas as pd
import pydynverse as pdv

# 测试用例1: 二分支结构数据


def get_test_wrap_data_bifurcation():
    # 测试样例数据, 与test_wrap_add_waypoints.py实际轨迹相同
    id = "test_wrap_add_branch_trajectory_bifurcation"
    cell_ids = ["a", "b", "c", "d", "e", "f"]
    branch_ids = ["A", "B", "C", "D"]
    branch_network = pd.DataFrame(
        columns=["from", "to"],
        data=[
            ["A", "B"],
            ["A", "C"],
            ["B", "D"],
        ],
    )
    branch_progressions = pd.DataFrame(
        columns=["cell_id", "branch_id", "percentage"],
        data=[
            ["a", "A", 0.0],
            ["b", "A", 0.8],
            ["c", "B", 0.2],
            ["d", "B", 1.0],
            ["e", "C", 0.2],  # 此处e细胞微调到分支C 0.2上
            ["f", "D", 0.2],
        ]
    )
    branches = pd.DataFrame(
        columns=["branch_id", "length", "directed"],
        data=[
            ["A", 1.0, True],
            ["B", 1.0, True],
            ["C", 1.0, True],
            ["D", 2.0, True],
        ]
    )

    # 封装数据
    dataset = pdv.wrap.wrap_data(
        id=id,
        cell_ids=cell_ids
    )
    return dataset, branch_network, branch_progressions, branches


# 测试用例1: 二分支结构执行测试
def test_wrap_add_branch_trajectory_bifurcation():

    dataset, branch_network, branch_progressions, branches = get_test_wrap_data_bifurcation()

    # 添加轨迹
    dataset = pdv.wrap.add_branch_trajectory(
        dataset=dataset,
        branch_ids=branches["branch_id"],
        branch_network=branch_network,
        branch_progressions=branch_progressions,
        branches=branches,
    )

    # 预期构造结果
    expected_milestone_network = pd.DataFrame(
        columns=["from", "to", "length", "directed"],
        data=[
            ["1", "2", 1.0, True],
            ["2", "3", 1.0, True],
            ["2", "4", 1.0, True],
            ["3", "5", 2.0, True],
        ]
    )
    expected_progressions = pd.DataFrame(
        columns=["cell_id", "from", "to", "percentage"],
        data=[
            ["a", "1", "2", 0.0],
            ["b", "1", "2", 0.8],
            ["c", "2", "3", 0.2],
            ["d", "2", "3", 1.0],
            ["e", "2", "4", 0.2],
            ["f", "3", "5", 0.2],
        ]
    )

    # 开始断言
    assert dataset["milestone_network"].equals(expected_milestone_network)
    assert dataset["progressions"].equals(expected_progressions)


# 测试用例2: 环结构数据
def get_test_wrap_data_loop():
    id = "test_wrap_add_branch_trajectory_loop"
    cell_ids = ["a", "b", "c", "d", "e", "f"]
    branch_ids = ["A", "B", "C"]
    branch_network = pd.DataFrame(
        columns=["from", "to"],
        data=[
            ["A", "B"],
            ["B", "C"],
            ["B", "A"],
            ["A", "C"],
        ]

    )
    branch_progressions = pd.DataFrame(
        columns=["cell_id", "branch_id", "percentage"],
        data=[
            ["a", "A", 0.0],
            ["b", "A", 0.5],
            ["c", "A", 1.0],
            ["d", "B", 0.1],
            ["e", "B", 0.2],
            ["f", "C", 0.3],
        ]
    )
    branches = pd.DataFrame(
        columns=["branch_id", "length", "directed"],
        data=[
            ["A", 1.0, True],
            ["B", 2.0, True],
            ["C", 3.0, True],
        ]
    )

    # 封装数据
    dataset = pdv.wrap.wrap_data(
        id=id,
        cell_ids=cell_ids
    )
    return dataset, branch_network, branch_progressions, branches


# 测试用例2: 环结构执行测试
def test_wrap_add_branch_trajectory_loop():
    dataset, branch_network, branch_progressions, branches = get_test_wrap_data_loop()

    # 添加轨迹
    dataset = pdv.wrap.add_branch_trajectory(
        dataset=dataset,
        branch_ids=branches["branch_id"],
        branch_network=branch_network,
        branch_progressions=branch_progressions,
        branches=branches,
    )

    # # 预期构造结果
    # expected_milestone_network = pd.DataFrame(
    #     columns=["from", "to", "length", "directed"],
    #     data=[
    #         ["1", "2", 1.0, True],
    #         ["2", "3", 1.0, True],
    #         ["2", "4", 1.0, True],
    #         ["3", "5", 2.0, True],
    #     ]
    # )
    # expected_progressions = pd.DataFrame(
    #     columns=["cell_id", "from", "to", "percentage"],
    #     data=[
    #         ["a", "1", "2", 0.0],
    #         ["b", "1", "2", 0.8],
    #         ["c", "2", "3", 0.2],
    #         ["d", "2", "3", 1.0],
    #         ["e", "2", "4", 0.2],
    #         ["f", "3", "5", 0.2],
    #     ]
    # )

    # # 开始断言
    # assert dataset["milestone_network"].equals(expected_milestone_network)
    # assert dataset["progressions"].equals(expected_progressions)


if __name__ == "__main__":
    pytest.main(["-v", __file__])
