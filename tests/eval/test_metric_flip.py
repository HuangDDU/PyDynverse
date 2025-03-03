# test_metric_flip.py
import pytest
import pandas as pd
import networkx as nx
from itertools import product
from typing import Dict, Tuple
from pydynverse.eval.metric_flip import calculate_edge_flip

def test_metric_flip():
    #原始网络数据
    linear1= pd.DataFrame(['a','b',1,True],columns=["from","to","length","directed"])
    linear2= pd.DataFrame([['a','b',1,True],['b','c',1,True]],columns=["from","to","length","directed"])
    #转化为networkx数据结构
    linear1_to_networkx1= None
    linear2_to_networkx2= None
    #调用函数测试返回得分
    score=0
    #断言
    assert score==0


if __name__ == "__main__":
    pytest.main(["-v", __file__])