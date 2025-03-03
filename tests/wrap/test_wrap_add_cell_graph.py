import pytest
import pydynverse as pdv

import pandas as pd
from ..test_util import compare_dataframes_closely


def test_wrap_add_cell_graph():
    dataset = None
    cell_ids = ["W", "X", "Y", "Z", "A", "WbX", "XcZ", "XeY", "ZfA", "a", "b", "c", "d", "e", "f"]
    dataset = pdv.wrap.wrap_data(
        id="test_wrap_add_cell_graph",
        cell_ids=cell_ids
    )
    cell_graph = pd.DataFrame(
        columns=["from", "to", "length", "directed"],
        data=[
            ["W", "WbX", 0.8, False],
            ["WbX", "X", 0.2, False],
            ["X", "XeY", 0.5, False],
            ["XeY", "Y", 0.5, False],
            ["X", "XcZ", 0.2, False],
            ["XcZ", "Z", 0.8, False],
            ["Z", "ZfA", 0.2, False],
            ["ZfA", "A", 0.8, False],
            ["W", "a", 0.5, False],
            ["WbX", "b", 0.5, False],
            ["XcZ", "c", 0.5, False],
            ["Z", "d", 0.5, False],
            ["XeY", "e", 0.2, False],
            ["ZfA", "f", 0.5, False],
        ]
    )
    cell_graph["directed"] = True  # 使用有向图，暂时比较简单

    to_keep = dict(
        W=True,
        X=True,
        Y=True,
        Z=True,
        A=True,
        WbX=True,
        XcZ=True,
        XeY=True,
        ZfA=True,
        a=False,
        b=False,
        c=False,
        d=False,
        e=False,
        f=False
    )
    to_keep = pd.Series(to_keep)

    # 执行
    trajectory = pdv.wrap.add_cell_graph(
        dataset=dataset,
        cell_graph=cell_graph,
        to_keep=to_keep,
        milestone_prefix="ML_",
    )

    # 预期输出
    expected_milestone_ids = [f"ML_{i}"for i in ["W", "X", "Y", "A"]]
    expected_milestone_network = pd.DataFrame(
        columns=["from", "to", "length", "directed"],
        data=[
            ["ML_W", "ML_X", 1, False],
            ["ML_X", "ML_Y", 1, False],
            ["ML_X", "ML_A", 2, False],
        ]
    )
    expected_milestone_network["directed"] = True  # 使用有向图，暂时比较简单

    expected_progressions = pd.DataFrame(
        columns=["cell_id", "from", "to", "percentage"],
        data=[
            ["W", "ML_W", "ML_X", 0],
            ["X", "ML_W", "ML_X", 1],
            ["Y", "ML_X", "ML_Y", 1],
            ["Z", "ML_X", "ML_A", 0.5],
            ["A", "ML_X", "ML_A", 1],
            ["WbX", "ML_W", "ML_X", 0.8],
            ["XcZ", "ML_X", "ML_A", 0.1],
            ["XeY", "ML_X", "ML_Y", 0.5],
            ["ZfA", "ML_X", "ML_A", 0.6],
            ["a", "ML_W", "ML_X", 0],
            ["b", "ML_W", "ML_X", 0.8],
            ["c", "ML_X", "ML_A", 0.1],
            ["d", "ML_X", "ML_A", 0.5],
            ["e", "ML_X", "ML_Y", 0.5],
            ["f", "ML_X", "ML_A", 0.6],
        ]
    )

    # assert
    assert trajectory["milestone_ids"] == expected_milestone_ids
    assert compare_dataframes_closely(trajectory["milestone_network"], expected_milestone_network, on_columns=["from", "to"])
    assert compare_dataframes_closely(trajectory["progressions"], expected_progressions, on_columns=["cell_id"])


if __name__ == "__main__":
    pytest.main(["-v", __file__])
