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
    
    # 构造里程碑网络
    milestone_network = pd.DataFrame(
        columns=["from", "to", "length", "directed"],
        data=[
            ["W", "X", 2, True],
            ["X", "Y", 3, True],
            ["X", "Z", 4, True],
        ]
    )
    
    # 构造分歧区域
    divergence_regions = pd.DataFrame(
        columns=["divergence_id", "milestone_id", "is_start"],
        data=[
            ["XYZ", "X", True],
            ["XYZ", "Y", False],
            ["XYZ", "Z", False]
        ]
    )
    
    # 构造里程碑百分比（长格式数据）
    milestone_percentages = pd.DataFrame(
        columns=["cell_id", "milestone_id", "percentage"],
        data=[
            ["a", "W", 0.9],
            ["a", "X", 0.1],
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
    
    # 构造表达矩阵，行索引为 cell_ids，列名固定为 ["gene_0", "gene_1", "gene_2"]
    expression = pd.DataFrame({
        "gene_0": [1, 2, 3, 4, 5],
        "gene_1": [2, 3, 4, 5, 6],
        "gene_2": [3, 4, 5, 6, 7],
    }, index=cell_ids)
    
    # 调用 wrap_data 时传入 feature_ids，确保表达矩阵列名与 feature_ids 一致
    dataset = pdv.wrap.wrap_data(id=id, cell_ids=cell_ids, feature_ids=list(expression.columns))
    # 如果 dataset 中没有 feature_info，则手动构造（索引使用 feature_ids）
    if "feature_info" not in dataset or dataset["feature_info"] is None:
        dataset["feature_info"] = pd.DataFrame(index=dataset["feature_ids"])
    
    # 调用 add_expression 将表达矩阵添加到数据集中
    # 此处 counts 也直接使用 expression（实际情况中 counts 可能与 expression 不同，但测试中可简化）
    dataset = add_expression(data=dataset, counts=expression, expression=expression)
    
    # 添加先验信息，用于 enrichment 测试，例如认为 "gene_0" 和 "gene_2" 重要
    dataset["prior_information"] = {"features_id": ["gene_0", "gene_2"]}
    
    # 整合其它数据
    test_data = {
        "dataset": dataset,
        "milestone_network": milestone_network,
        "divergence_regions": divergence_regions,
        "milestone_percentages": milestone_percentages,
    }
    
    # 取出各项数据
    dataset = test_data["dataset"]
    milestone_network = test_data["milestone_network"]
    divergence_regions = test_data["divergence_regions"]
    milestone_percentages = test_data["milestone_percentages"]
    
    # 调用 add_trajectory 将轨迹信息添加到数据集中
    test_trajectory = pdv.wrap.add_trajectory(
        dataset,
        milestone_network=milestone_network,
        divergence_regions=divergence_regions,
        milestone_percentages=milestone_percentages
    )
    return test_trajectory

def test_function():
    # 构造数据集并检查生成的轨迹对象是否包含必要信息
    trajectory = get_data()
    # 检查轨迹对象中包含 expression、milestone_percentages、milestone_ids、cell_ids 及 prior_information
    assert "expression" in trajectory, "轨迹对象缺少 expression"
    assert "milestone_percentages" in trajectory, "轨迹对象缺少 milestone_percentages"
    assert "milestone_ids" in trajectory, "轨迹对象缺少 milestone_ids"
    assert "cell_ids" in trajectory, "轨迹对象缺少 cell_ids"
    assert "prior_information" in trajectory, "轨迹对象缺少 prior_information"
    return trajectory

def test_metric_featureimp():
    # 使用真实接口进行整体逻辑测试
    trajectory = test_function()
    
    # 调用 calculate_overall_feature_importance
    overall_imp = calculate_overall_feature_importance(
        trajectory=trajectory,
        expression_source="expression",
        fi_method=fi_ranger_rf_lite(),
        verbose=False
    )
    # 检查返回结果为 DataFrame 且包含 'feature_id' 和 'importance'
    assert isinstance(overall_imp, pd.DataFrame)
    expected_cols = {"feature_id", "importance"}
    assert set(overall_imp.columns) == expected_cols, f"预期列 {expected_cols}, 得到 {set(overall_imp.columns)}"
    
    # 调用 calculate_featureimp_cor（这里使用相同轨迹作为 dataset 和 prediction）
    cor_result = calculate_featureimp_cor(
        dataset=trajectory,
        prediction=trajectory,
        expression_source="expression",
        fi_method=fi_ranger_rf_lite()
    )
    # 检查返回的字典中包含预期键，且数值在 0 到 1 之间
    assert isinstance(cor_result, dict)
    assert "featureimp_cor" in cor_result and "featureimp_wcor" in cor_result
    assert 0 <= cor_result["featureimp_cor"] <= 1
    assert 0 <= cor_result["featureimp_wcor"] <= 1

    # 调用 calculate_featureimp_enrichment
    enrich_result = calculate_featureimp_enrichment(
        dataset=trajectory,
        prediction=trajectory,
        expression_source="expression",
        fi_method=fi_ranger_rf_lite()
    )
    # 检查 enrichment 返回结果格式及数值范围
    assert isinstance(enrich_result, dict)
    assert "featureimp_ks" in enrich_result and "featureimp_wilcox" in enrich_result
    assert 0 <= enrich_result["featureimp_ks"] <= 1
    assert 0 <= enrich_result["featureimp_wilcox"] <= 1

if __name__ == "__main__":
    pytest.main(["-v", __file__])
