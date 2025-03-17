import pytest
import pydynverse as pdv

import pandas as pd


def test_calculate_metrics1():
    linear = pd.DataFrame(
        columns=["from", "to", "length", "directed"],
        data=[
            ["A", "B", 1, True,],
            ["B", "C", 2, True,],
            ["C", "D", 3, True,],
        ],
    ) # 会对线性简化
    bifurcating = pd.DataFrame(
        data=[
            ["A", "B", 1, True,],
            ["B", "C", 2, True,],
            ["B", "D", 3, True,],
        ],
        columns=["from", "to", "length", "directed"]
    )
    unsimplified_score = pdv.eval.calculate_edge_flip(linear, bifurcating, simplify=False)
    # simplified_score = pdv.eval.calculate_edge_flip(linear, bifurcating, simplify=True)

    expected_unsimplified_score = 1 - 2/4
    # expected_simplified_score = 0

    assert unsimplified_score == expected_unsimplified_score
    # assert simplified_score == expected_simplified_score

def test_calculate_metrics2():
    # TODO: 
    pass

if __name__ == "__main__":
    pytest.main(["-v", __file__])
