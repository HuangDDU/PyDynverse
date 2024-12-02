import pytest

import pandas as pd
from pydynverse.wrap import wrap_data, is_data_wrapper


def test_wrap_data():
    # 最常见的情况，给定cell_ids, cell_info（adata.obs中细胞编号与批次信息或其他注释），feature_ids（adata.var中的基因名）
    id = "a"
    cell_ids = ["a", "b", "c", "d", "e"]
    cell_info = pd.DataFrame({
        "cell_ids": cell_ids,
        "info1": [f"info1_{cell_id}" for cell_id in cell_ids],
        "info2": [f"info2_{cell_id}" for cell_id in cell_ids],
    })
    feature_ids = ["g1", "g2", "g3"]
    expected_feature_info = pd.DataFrame(index= feature_ids) # 预期输出的feature_info DataFrame

    dataset = wrap_data(id=id, cell_ids=cell_ids, cell_info=cell_info, feature_ids=feature_ids)

    assert is_data_wrapper(dataset), "wrapper data failed"
    
    assert dataset["cell_ids"]==cell_ids, "cell_ids init failed"
    assert dataset["cell_info"].equals(cell_info), "cell_info init failed"
    assert dataset["feature_ids"]==feature_ids, "feature_ids init failed"
    assert dataset["feature_info"].equals(expected_feature_info), "feature_info init failed"



if __name__ == "__main__":
    pytest.main()
