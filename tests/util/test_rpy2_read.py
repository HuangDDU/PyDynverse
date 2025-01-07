import pytest
import pandas as pd

import rpy2.robjects as ro
from pydynverse.util import rpy2_read # 导入对应的装饰器, 数据自动转换

from scipy.sparse import csc_matrix

def test_rpy2_read():
    # 读取quickstart示例数据
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


def test_rpy2_read_sim():
    # 读取模拟数据
    r_script = """
    dataset <- readRDS("/home/huang/RCode/scrna_tools/dynbenchmark/data/synthetic/dyntoy/bifurcating_1.rds") # 二分支结构
    dataset
    """
    # 这里数据读取有点错乱
    dataset = ro.r(r_script) # 读取的为dcgMatrix，读取的为Matrix
    assert isinstance(dataset["counts"], csc_matrix), "dataset['counts'] is not a csc_matrix"
    assert isinstance(dataset["expression"], csc_matrix), "dataset['expression'] is not a csc_matrix"

if __name__ == "__main__":
    pytest.main(["-v", __file__])
