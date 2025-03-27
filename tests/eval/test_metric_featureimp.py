import numpy as np
import pandas as pd
import pydynverse as pdv
import pytest
from pydynverse.feature.fi_methods import fi_ranger_rf_lite
from pydynverse.feature.calculate_overall_feature_importance import calculate_overall_feature_importance
from pydynverse.eval.metric_featureimp import (
    calculate_featureimp_cor,
    calculate_featureimp_enrichment,
)

from pydynverse.wrap.wrap_add_expression import add_expression

def get_data():
    id = "test_metric_featureimp"
    cell_ids = ["a", "b", "c", "d", "e"]
    milestone_ids = ["W", "X", "Y", "Z"]
    milestone_network = pd.DataFrame(
        columns=["from", "to", "length", "directed"],
        data=[
            ["W", "X", 2, True],
            ["X", "Y", 3, True],
            ["X", "Z", 4, True],
        ]
    )
    divergence_regions = pd.DataFrame(
        columns=["divergence_id", "milestone_id", "is_start"],
        data=[
            ["XYZ", "X", True],
            ["XYZ", "Y", False],
            ["XYZ", "Z", False]
        ]
    )
    milestone_percentages = pd.DataFrame(
        columns=["cell_id", "milestone_id", "percentage"],
        data=[
            ["a", "W", 0.9],
            ["a","X",0.1],
            ["b", "W", 0.2],
            ["b", "X", 0.8],
            ["c", "X", 0.8],
            ["c", "Z", 0.2],
            ["d", "Z", 0.1],
            ["d", "X", 0.2],
            ["d", "Y", 0.7],
            ["e", "X", 0.3],
            ["e", "Y", 0.2],
            ["e", "Z", 0.5],
        ]
    )
    dataset = pdv.wrap.wrap_data(id=id, cell_ids=cell_ids)
    #构造expression
    expression=None
    dataset=add_expression(data=dataset,counts=expression,expression = expression)


    test_data = {
        "dataset": dataset,
        "milestone_network": milestone_network,
        "divergence_regions": divergence_regions,
        "milestone_percentages": milestone_percentages,
    }


    dataset=test_data["dataset"]
    milestone_network=test_data["milestone_network"]
    divergence_regions=test_data["divergence_regions"]
    milestone_percentages=test_data["milestone_percentages"]
    # 添加轨迹
    test_trajectory = pdv.wrap.add_trajectory(
        dataset,
        milestone_network=milestone_network,
        divergence_regions=divergence_regions,
        milestone_percentages=milestone_percentages
    )
    return test_trajectory

def test_function():
    #构建数据集
    assert 1==1
    #函数功能测试
    

def test_metric_featureimp():
    #最终逻辑测试
    assert 1==1


if __name__ == "__main__":
    pytest.main(["-v", __file__])