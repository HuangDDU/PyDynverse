import pytest
import pydynverse as pdv
import pandas as pd


def get_test_wrap_data():
    # 输入数据
    id = "test_add_end_state_probabilities"
    cell_ids = ["a", "aa", "b", "bb", "c", "cc"]
    end_state_ids = ["A", "B", "C"]
    end_state_probabilities = pd.DataFrame(
        columns=["cell_id", "A", "B", "C"],
        data=[
            ["a", .5, 0, 0],
            ["aa", 1, 0, 0],
            ["b", 0, .5, 0],
            ["bb", 0, 1, 0],
            ["c", 0, 0, .5],
            ["cc", 0, 0, 1],
        ]
    )
    pseudotime = [.5, 1, .5, 1, .5, 1]
    pseudotime = pd.Series(pseudotime, index=cell_ids)
    dataset = pdv.wrap.wrap_data(cell_ids=cell_ids, id=id)

    test_wrap_data = {
        "dataset": dataset,
        "end_state_probabilities": end_state_probabilities,
        "end_state_ids": end_state_ids,
        "pseudotime": pseudotime,

    }

    return test_wrap_data


# 测试样例，3个状态构成星型结构
def test_add_end_state_probabilities_3_states():
    # 输入数据
    test_wrap_data = get_test_wrap_data()
    dataset = test_wrap_data["dataset"]
    end_state_probabilities = test_wrap_data["end_state_probabilities"]
    end_state_ids = test_wrap_data["end_state_ids"]
    pseudotime = test_wrap_data["pseudotime"]

    # 执行
    trajectory = pdv.wrap.add_end_state_probabilities(
        dataset=dataset,
        end_state_probabilities=end_state_probabilities,
        pseudotime=pseudotime,
    )

    # 预期输出
    start_milestone_id = "milestone_begin"
    milestone_ids = [start_milestone_id] + end_state_ids
    expected_milestone_network = pd.DataFrame({
        "from": start_milestone_id,
        "to": end_state_ids,
        "length": 1,
        "directed": True
    })
    expected_divergence_regions = pd.DataFrame({
        "milestone_id": milestone_ids,
        "divergence_id": "D",
        "is_start": pd.Series(milestone_ids) == start_milestone_id
    })
    scaled_pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())
    expected_progressions = end_state_probabilities.melt(id_vars=["cell_id"], var_name="to", value_name="percentage")
    expected_progressions["from"] = start_milestone_id
    expected_progressions["percentage"] = expected_progressions.groupby("cell_id")["percentage"].transform(lambda x: x / x.sum() * scaled_pseudotime[x.name])  # 缩放使其之和为1，暂时不理解这个
    expected_progressions = expected_progressions[["cell_id", "from", "to", "percentage"]]

    assert trajectory["milestone_network"].equals(expected_milestone_network)
    assert trajectory["divergence_regions"].equals(expected_divergence_regions)
    assert trajectory["progressions"].equals(expected_progressions)


# 不指定状态，就是线性轨迹，直接使用伪时间
def test_add_end_state_probabilities_without_state():
    # 输入数据
    test_wrap_data = get_test_wrap_data()
    dataset = test_wrap_data["dataset"]
    end_state_probabilities = test_wrap_data["end_state_probabilities"]
    pseudotime = test_wrap_data["pseudotime"]
    end_state_probabilities = end_state_probabilities["cell_id"].to_frame() # 没有终端状态

    
    trajectory = pdv.wrap.add_end_state_probabilities(
        dataset=dataset,
        end_state_probabilities=end_state_probabilities,
        pseudotime=pseudotime,
    )

    # 预期输出，相当于直接调用线性轨迹
    expected_trajectory = pdv.wrap.add_linear_trajectory(
        dataset=dataset,
        pseudotime=pseudotime,
        directed=True,
    )

    assert trajectory["milestone_network"].equals(expected_trajectory["milestone_network"])
    assert trajectory["divergence_regions"].equals(expected_trajectory["divergence_regions"])
    assert trajectory["progressions"].equals(expected_trajectory["progressions"])
    


if __name__ == "__main__":
    pytest.main(["-v", __file__])
