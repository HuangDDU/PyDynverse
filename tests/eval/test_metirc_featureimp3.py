import numpy as np
import pandas as pd
import pytest
from pydynverse.feature.fi_methods import fi_ranger_rf_lite
from pydynverse.feature.calculate_overall_feature_importance import calculate_overall_feature_importance
from pydynverse.eval.metric_featureimp import (
    calculate_featureimp_cor,
    calculate_featureimp_enrichment,
)

# 构造一个固定数据集，保证数值可重复
@pytest.fixture
def dummy_trajectory_dataset():
    """
    构造一个 dataset 轨迹，包含：
       cell_ids: 固定的细胞列表
       expression: 固定的表达矩阵（4个细胞 x 5个基因），数值固定
       milestone_percentages: 每个细胞归属到一个里程碑（简单构造）
       milestone_ids: ["m1", "m2"]
       prior_information: 包含先验特征列表（用于 enrichment）
    """
    cell_ids = ["cell1", "cell2", "cell3", "cell4"]
    # 构造表达矩阵，数据固定
    expr = pd.DataFrame({
        "gene0": [1, 2, 3, 4],
        "gene1": [2, 2, 2, 2],
        "gene2": [1, 3, 5, 7],
        "gene3": [4, 5, 6, 7],
        "gene4": [1, 1, 1, 1]
    }, index=cell_ids)
    # 每个细胞只属于一个里程碑（例如前两细胞 m1，后两细胞 m2）
    milestone_percentages = pd.DataFrame({
        "cell_id": cell_ids,
        "milestone_id": ["m1", "m1", "m2", "m2"],
        "percentage": [1.0, 1.0, 1.0, 1.0]
    })
    return {
        "cell_ids": cell_ids,
        "expression": expr,
        "milestone_percentages": milestone_percentages,
        "milestone_ids": ["m1", "m2"],
        "pydynwrap:with_trajectory": True,
        "prior_information": {"features_id": ["gene0", "gene2", "gene4"]}  # 假定先验认为这3个基因重要
    }

@pytest.fixture
def dummy_trajectory_prediction():
    """
    构造一个 prediction 轨迹，与 dataset 类似，但对表达矩阵做线性变换（例如乘以2），
    使得两者在整体特征重要性上存在线性关系。
    """
    cell_ids = ["cell1", "cell2", "cell3", "cell4"]
    expr = pd.DataFrame({
        "gene0": [2, 4, 6, 8],
        "gene1": [4, 4, 4, 4],
        "gene2": [2, 6, 10, 14],
        "gene3": [8, 10, 12, 14],
        "gene4": [2, 2, 2, 2]
    }, index=cell_ids)
    milestone_percentages = pd.DataFrame({
        "cell_id": cell_ids,
        "milestone_id": ["m1", "m1", "m2", "m2"],
        "percentage": [1.0, 1.0, 1.0, 1.0]
    })
    return {
        "cell_ids": cell_ids,
        "expression": expr,
        "milestone_percentages": milestone_percentages,
        "milestone_ids": ["m1", "m2"],
        "pydynwrap:with_trajectory": True
    }

# 测试 calculate_overall_feature_importance
def test_calculate_overall_feature_importance_integration(dummy_trajectory_dataset):
    # 使用默认的 fi_ranger_rf_lite
    overall_imp = calculate_overall_feature_importance(
        trajectory=dummy_trajectory_dataset,
        expression_source="expression",
        fi_method=fi_ranger_rf_lite(),
        verbose=False
    )
    # 检查返回结果为 DataFrame 且包含 feature_id 和 importance 两列
    assert isinstance(overall_imp, pd.DataFrame)
    expected_cols = {"feature_id", "importance"}
    assert set(overall_imp.columns) == expected_cols, f"Expected columns {expected_cols}, got {set(overall_imp.columns)}"
    # 检查结果非空且 importance 值均为非负数
    assert not overall_imp.empty
    assert (overall_imp["importance"] >= 0).all()
    # 检查排序：第一行的 importance 应为最大值
    sorted_importance = overall_imp["importance"].sort_values(ascending=False).values
    np.testing.assert_allclose(overall_imp["importance"].values, sorted_importance, rtol=1e-3)


# 测试 calculate_featureimp_cor（整体函数）

def test_calculate_featureimp_cor_integration(dummy_trajectory_dataset, dummy_trajectory_prediction):
    result = calculate_featureimp_cor(
        dataset=dummy_trajectory_dataset,
        prediction=dummy_trajectory_prediction,
        expression_source="expression",
        fi_method=fi_ranger_rf_lite()
    )
    # 检查返回字典中包含 featureimp_cor 和 featureimp_wcor，且数值在合理范围内（0到1之间）
    assert "featureimp_cor" in result and "featureimp_wcor" in result
    assert 0 <= result["featureimp_cor"] <= 1
    assert 0 <= result["featureimp_wcor"] <= 1


# 测试 calculate_featureimp_enrichment（整体函数）

def test_calculate_featureimp_enrichment_integration(dummy_trajectory_dataset, dummy_trajectory_prediction):
    result = calculate_featureimp_enrichment(
        dataset=dummy_trajectory_dataset,
        prediction=dummy_trajectory_prediction,
        expression_source="expression",
        fi_method=fi_ranger_rf_lite()
    )
    # 如果 prediction 中 cell_id 数量满足条件，返回应为字典，且 ks 与 wilcox 值在 [0,1] 内
    assert "featureimp_ks" in result and "featureimp_wilcox" in result
    assert 0 <= result["featureimp_ks"] <= 1
    assert 0 <= result["featureimp_wilcox"] <= 1


# 测试 prediction 不满足条件时（例如 cell_id 数量不足）返回默认值

def test_featureimp_invalid_prediction(dummy_trajectory_dataset):
    invalid_prediction = {
        "milestone_percentages": pd.DataFrame({
            "cell_id": ["c1", "c1"],  # 唯一 cell_id 数量不足3
            "milestone_id": ["m1", "m1"],
            "percentage": [1.0, 1.0]
        })
    }
    result_cor = calculate_featureimp_cor(dummy_trajectory_dataset, invalid_prediction)
    result_enrich = calculate_featureimp_enrichment(dummy_trajectory_dataset, invalid_prediction)
    assert result_cor == {'featureimp_cor': 0, 'featureimp_wcor': 0}
    assert result_enrich == {'featureimp_ks': 0, 'featureimp_wilcox': 0}

if __name__ == "__main__":
    pytest.main(["-v", __file__])
