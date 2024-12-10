import pytest
import pandas as pd
import pydynverse as pdv


def test_wrap_add_branch_trajectory():
    # 测试样例数据, 与test_wrap_add_waypoints.py实际轨迹相同
    id = "test_wrap_add_branch_trajectory"
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
    wr_orig = pdv.wrap.wrap_data(
        id=id,
        cell_ids=cell_ids
    )

    # 添加轨迹
    wr = pdv.wrap.add_branch_trajectory(
        dataset=wr_orig,
        branch_ids=branch_ids,
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
    assert wr["milestone_network"].equals(expected_milestone_network)
    assert wr["progressions"].equals(expected_progressions)


if __name__ == "__main__":
    pytest.main(["-v", __file__])
