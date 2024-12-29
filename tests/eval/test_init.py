import pytest
import pydynverse as pdv

import pandas as pd


def test_init():
    metrics = pdv.eval.metrics
    assert isinstance(metrics, pd.DataFrame)
    assert set(["edge_flip", "isomorphic"]) <= set(metrics["metric_id"])


if __name__ == "__main__":
    pytest.main(["-v", __file__])
