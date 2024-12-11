import pytest

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pydynverse as pdv
# from ..wrap.test_wrap_add_waypoints import get_test_wrap_data # 这种相对路径导入失败
from tests.wrap.test_wrap_add_waypoints import get_test_wrap_data

pdv.logger.setLevel("DEBUG")

def test_plot_dimred():
    # 测试样例数据
    dataset, milestone_network, divergence_regions, milestone_percentages = get_test_wrap_data()

    # 重新封装成带着AnnData的数据
    cell_ids = dataset["cell_ids"]
    feature_ids = ["G1", "G2"]
    expression = np.zeros((len(cell_ids), len(feature_ids)))
    counts = expression.copy()
    dataset = pdv.wrap.wrap_expression(
        expression=expression,
        counts=counts,
        cell_ids=cell_ids,
        feature_ids=feature_ids
    )
    trajectory = pdv.wrap.add_trajectory(
        dataset,
        milestone_network=milestone_network,
        divergence_regions=divergence_regions,
        milestone_percentages=milestone_percentages
    )

    # 绘图需要的额外数据
    grouping = [1, 1, 2, 2, 2, 3]
    space = pd.DataFrame(
        columns=["comp_1", "comp_2"],
        index=dataset["cell_ids"],
        data=[
            [0, 10],
            [8, 10],
            [12, 12],
            [20, 20],
            [15, 16],
            [22, 20],
        ]
    )
    space["cell_id"] = space.index

    # 添加轨迹
    trajectory = pdv.wrap.add_trajectory(
        dataset,
        milestone_network=milestone_network,
        divergence_regions=divergence_regions,
        milestone_percentages=milestone_percentages
    )

    # 绘图
    pdv.plot.plot_dimred(trajectory, grouping=grouping, dimred=space)
    plt.savefig("img/test_plot_dimred.png")


if __name__ == "__main__":
    pytest.main(["-v", __file__])