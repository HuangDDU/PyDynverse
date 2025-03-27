import numpy as np
import pandas as pd
import pytest
from pydynverse.feature.calculate_overall_feature_importance import calculate_overall_feature_importance

# 定义一个 dummy milestone 特征重要性函数，用于模拟 calculate_milestone_feature_importance 的返回结果
def dummy_milestone_feature_importance(trajectory, expression_source, fi_method, verbose):
    # 模拟返回数据：假设有两个里程碑（"m1" 和 "m2"），每个里程碑计算出两个特征的 importance
    data = [
        {"milestone_id": "m1", "feature_id": "gene0", "importance": 0.8},
        {"milestone_id": "m1", "feature_id": "gene1", "importance": 0.6},
        {"milestone_id": "m2", "feature_id": "gene0", "importance": 0.7},
        {"milestone_id": "m2", "feature_id": "gene1", "importance": 0.9},
    ]
    return pd.DataFrame(data)

# 构造一个 dummy 轨迹数据（trajectory）
@pytest.fixture
def dummy_trajectory():
    # 细胞 ID 列表（虽然此处不会实际使用到表达数据，因为我们将通过 monkeypatch 替换 milestone 计算）
    cell_ids = ["cell1", "cell2", "cell3"]
    # 构造一个简单的表达矩阵（行为细胞，列为基因）
    expr = pd.DataFrame(
        np.random.rand(3, 2),
        index=cell_ids,
        columns=["gene0", "gene1"]
    )
    # 构造里程碑百分比的长格式数据（内容不会真正参与计算）
    milestone_percentages = pd.DataFrame({
        "cell_id": cell_ids * 2,
        "milestone_id": ["m1", "m1", "m1", "m2", "m2", "m2"],
        "percentage": [0.7, 0.5, 0.9, 0.3, 0.5, 0.1]
    })
    return {
        "cell_ids": cell_ids,
        "expression": expr,
        "milestone_percentages": milestone_percentages,
        "milestone_ids": ["m1", "m2"],
        "pydynwrap:with_trajectory": True
    }

def test_calculate_overall_feature_importance(monkeypatch, dummy_trajectory):
    # 使用 monkeypatch 替换 calculate_milestone_feature_importance，
    # 确保当 calculate_overall_feature_importance 内部调用它时，返回我们预定义的 dummy 数据
    monkeypatch.setattr(
        "pydynverse.feature.calculate_overall_feature_importance.calculate_milestone_feature_importance",
        dummy_milestone_feature_importance
    )

    # 调用 calculate_overall_feature_importance
    # 此处我们传入 dummy_fi_method 的占位函数（实际不会用到，因为 milestone 计算被替换掉了）
    dummy_fi_method = {"fun": lambda X, y, verbose=False: {}}
    result = calculate_overall_feature_importance(
        trajectory=dummy_trajectory,
        expression_source="expression",
        fi_method=dummy_fi_method,
        verbose=False
    )

    # 根据 dummy_milestone_feature_importance 的返回数据，
    # gene0 的平均 importance = (0.8 + 0.7) / 2 = 0.75
    # gene1 的平均 importance = (0.6 + 0.9) / 2 = 0.75
    expected = pd.DataFrame({
        "feature_id": ["gene0", "gene1"],
        "importance": [0.75, 0.75]
    })

    # 对结果按照 feature_id 排序后再比较
    result_sorted = result.sort_values("feature_id").reset_index(drop=True)
    expected_sorted = expected.sort_values("feature_id").reset_index(drop=True)
    pd.testing.assert_frame_equal(result_sorted, expected_sorted)

if __name__ == "__main__":
    pytest.main(["-v", __file__])
