import pytest

import pandas as pd
from pydynverse.wrap import wrap_expression, is_wrapper_with_expression


def test_wrap_add_expression():
    id = "a"
    cell_ids = ["a", "b", "c", "d", "e"]
    cell_info = pd.DataFrame({
        "cell_ids": cell_ids,
        "info1": [f"info1_{cell_id}" for cell_id in cell_ids],
        "info2": [f"info2_{cell_id}" for cell_id in cell_ids],
    })
    feature_ids = ["g1", "g2", "g3"]

    expression = [
        []
    ]

    wrap_expression()
    dataset = wrap_data(id=id, cell_ids=cell_ids, cell_info=cell_info, feature_ids=feature_ids)

    assert is_data_wrapper(dataset), "wrapper data failed"
    
    assert dataset["cell_ids"]==cell_ids, "cell_ids init failed"
    assert dataset["cell_info"].equals(cell_info), "cell_info init failed"
    assert dataset["feature_ids"]==feature_ids, "feature_ids init failed"
    assert dataset["feature_info"].equals(expected_feature_info), "feature_info init failed"



if __name__ == "__main__":
    pytest.main()