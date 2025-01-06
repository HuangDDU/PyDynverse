import pytest
import pydynverse as pdv

import pandas as pd


def test_calculate_metrics():
    # TODO: 构造数据集获得整体评分
    id = "directed_linear"
    cell_ids = ["a", "b", "c", "d", "e"]
    dataset_raw = pdv.wrap.wrap_data(id=id, cell_ids=cell_ids)  # 只有id

    # 参考milestone
    milestone_ids = ["A", "B", "C", "D"]
    milestone_network = pd.DataFrame(
        data=[
            ["A", "B", 1, True],
            ["B", "C", 1, True],
            ["C", "D", 1, True]
        ],
        columns=["from", "to", "length", "directed"],
    )
    progressions = pd.DataFrame(
        data=[
            ["a", "A", "B", 0.3],
            ["b", "A", "B", 0.6],
            ["c", "B", "C", 0.2],
            ["d", "B", "C", 0.8],
            ["e", "C", "D", 0.4],
        ],
        columns=["cell_id", "from", "to", "percentage"]
    )
    dataset = pdv.wrap.add_trajectory(
        dataset=dataset_raw.copy(),
        milestone_ids=milestone_ids,
        milestone_network=milestone_network,
        progressions=progressions
    )
    # 模型预测出的milestone
    milestone_ids = ["A", "B", "C", "D"]
    milestone_network = pd.DataFrame(
        data=[
            ["A", "B", 1, True],
            ["B", "C", 1, True],
            ["B", "D", 1, True]
        ],
        columns=["from", "to", "length", "directed"],
    )
    progressions = pd.DataFrame(
        data=[
            ["a", "A", "B", 0.3],
            ["b", "A", "B", 0.7],
            ["c", "B", "C", 0.2],
            ["d", "B", "C", 0.8],
            ["e", "B", "D", 0.5],
        ],
        columns=["cell_id", "from", "to", "percentage"]
    )
    model = pdv.wrap.add_trajectory(
        dataset=dataset_raw.copy(),
        milestone_ids=milestone_ids,
        milestone_network=milestone_network,
        progressions=progressions
    )

    summary_dict_self = pdv.eval.calculate_metrics(dataset, dataset)  # 自身计算指标，一定是指标最大值
    summary_dict = pdv.eval.calculate_metrics(dataset, model)
    assert summary_dict_self["isomorphic"] == 1
    assert summary_dict["isomorphic"] == 0


if __name__ == "__main__":
    pytest.main(["-v", __file__])
