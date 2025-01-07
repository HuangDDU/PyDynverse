import pytest
import pydynverse as pdv

import pandas as pd

from ..wrap.test_wrap_add_waypoints import get_test_wrap_data


def test_mix_colors():
    dataset, milestone_network, divergence_regions, milestone_percentages = get_test_wrap_data()  # 复用test_wrap_add_waypoints数据
    milestone_colors = pd.DataFrame(
        data=[
            [1, 0, 0],
            [1, 0.5, 0],
            [0, 0, 1],
            [1, 1, 0],
            [0, 1, 0],
        ],
        index=["W", "X", "Y", "Z", "A"],
        columns=["R", "G", "B"],  # 列名在脚本中也没有给出来
    )

    # 测试单个细胞b
    sc_milestone_percentages = milestone_percentages[milestone_percentages["cell_id"] == "b"]
    sc_color = pdv.plot.mix_colors(sc_milestone_percentages["milestone_id"], sc_milestone_percentages["percentage"], milestone_colors)
    assert sc_color == "#ff6600"

    # 测试多个细胞
    cell_colors = milestone_percentages.groupby("cell_id").apply(lambda x: pdv.plot.mix_colors(x["milestone_id"], x["percentage"], milestone_colors))
    assert cell_colors.iloc[1] == "#ff6600"

if __name__ == "__main__":
    pytest.main(["-v", __file__])
