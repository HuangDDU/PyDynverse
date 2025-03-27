import pytest
import pydynverse as pdv

import numpy as np
import pandas as pd
from .test_wrap_add_waypoints import get_test_wrap_data as get_test_wrap_data_ref


def get_test_wrap_data():
    test_wrap_data = get_test_wrap_data_ref()
    dataset = test_wrap_data["dataset"]
    cell_ids = dataset["cell_ids"]
    feature_ids = dataset["feature_ids"]
    feature_ids = feature_ids if feature_ids is not None else ["g1", "g2", "g3"]
    expression = np.zeros((len(cell_ids), len(feature_ids)))  # 此处需要添加expression
    counts = expression.copy()  # count与expression相同
    test_wrap_data = {
        "cell_ids": cell_ids,
        "feature_ids": feature_ids,
        "cell_info": dataset["cell_info"],
        "feature_info": dataset["feature_info"],
        "counts": counts,
        "expression": expression,
        "dataset": dataset,
    }
    return test_wrap_data


def test_add_expression():
    # 对已有的dataset添加表达矩阵
    test_wrap_data = get_test_wrap_data()
    dataset = test_wrap_data["dataset"]
    expression = test_wrap_data["expression"]
    counts = test_wrap_data["counts"]
    pdv.wrap.add_expression(dataset, counts=counts, expression=expression)

    assert "counts" in dataset
    assert "expression" in dataset


def test_wrap_expression():
    # 创建dataset并添加表达矩阵
    test_wrap_data = get_test_wrap_data()
    dataset = test_wrap_data["dataset"]
    expression = test_wrap_data["expression"]
    counts = test_wrap_data["counts"]
    cell_ids = test_wrap_data["cell_ids"]
    feature_ids = test_wrap_data["feature_ids"]

    dataset = pdv.wrap.wrap_expression(
        expression=expression,
        counts=counts,
        cell_ids=cell_ids,
        feature_ids=feature_ids
    )

    assert pdv.wrap.is_wrapper_with_expression(dataset)

if __name__ == "__main__":
    pytest.main()
