import pytest
import pandas as pd

import rpy2.robjects as ro
from pydynverse.util import rpy2_read # 导入对应的装饰器

def test_rpy2_read():
    # 测试从R语言到Python的自动转换
    r_script = """
    library(dyno)
    data("fibroblast_reprogramming_treutlein")
    fibroblast_reprogramming_treutlein
    """
    dataset = ro.r(r_script)
    assert type(dataset)==dict, "dataset read object is not a dict"
    assert type(dataset["cell_ids"])==list, "dataset's attibute 'cell_ids' is not a list"
    assert type(dataset["feature_info"])==pd.DataFrame, "dataset's attibute 'feature_info' is not a data frame"
    assert type(dataset["grouping"])==list, "dataset's attibute 'grouping' is not a list"
    # count和expression为字典，且其键为csc, cell_ids, feature_id
    assert type(dataset["counts"])==dict, "dataset's attibute 'count' is not a dict"
    assert type(dataset["expression"])==dict, "dataset's attibute 'expression' is not a dict"

if __name__ == "__main__":
    pytest.main(["-v", __file__])
