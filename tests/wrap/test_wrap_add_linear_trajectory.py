import pandas as pd
import pytest
import pydynverse as pdv


def get_test_wrap_data():
    id = "test_wrap_add_linear_trajectory"
    cell_ids = ["a", "b", "c", "d", "e", "f"]
    pseudotime = [0.0, 0.1, 0.4, 0.5, 0.8, 1.0]
    dataset = pdv.wrap.wrap_data(cell_ids=cell_ids, id=id)
    pdv.wrap.add_linear_trajectory(dataset, pseudotime)
    return dataset, cell_ids, pseudotime


def test_wrap_add_linear_trajectory():
    dataset, cell_ids, pseudotime = get_test_wrap_data()

    expected_milestone_ids = ["milestone_begin", "milestone_end"]
    expected_milestone_network = pd.DataFrame({
        "from": "milestone_begin",
        "to": "milestone_end",
        "length": 1,
        "directed": False,
    }, index=[0])
    expected_progressions = pd.DataFrame({
        "cell_id": cell_ids,
        "from": "milestone_begin",
        "to": "milestone_end",
        "percentage": pseudotime,
    })

    # 构造的milestone_network和progressions与预期对比
    assert dataset["milestone_ids"] == expected_milestone_ids
    assert dataset["milestone_network"].equals(expected_milestone_network)
    assert dataset["progressions"].equals(expected_progressions)

    assert dataset["directed"]  # 暂时都默认为True


if __name__ == "__main__":
    pytest.main(["-v", __file__])
