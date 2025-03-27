import numpy as np
import pandas as pd
import pytest
from pydynverse.feature.calculate_overall_feature_importance import calculate_overall_feature_importance

def dummy_milestone_feature_importance_complex(trajectory, expression_source, fi_method, verbose):
    # 获取所有里程碑和特征列表
    milestone_ids = trajectory["milestone_ids"]  # ["m1", "m2", "m3"]
    feature_ids = ["gene0", "gene1", "gene2", "gene3", "gene4"]
    
    # 生成全量组合（所有特征 x 所有里程碑），默认重要性为0
    data = []
    for milestone in milestone_ids:
        for feature in feature_ids:
            data.append({
                "milestone_id": milestone,
                "feature_id": feature,
                "importance": 0.0  # 默认补零
            })
    
    # 更新实际有值的条目
    updates = [
        ("m1", "gene0", 0.9),
        ("m1", "gene1", 0.6),
        ("m1", "gene2", 0.0),
        ("m2", "gene0", 0.7),
        ("m2", "gene1", 0.8),
        ("m2", "gene3", 0.4),
        ("m3", "gene0", 0.5),
        ("m3", "gene2", 0.3),
        ("m3", "gene4", 0.6),
    ]
    for milestone, feature, imp in updates:
        for entry in data:
            if entry["milestone_id"] == milestone and entry["feature_id"] == feature:
                entry["importance"] = imp
                break
    
    return pd.DataFrame(data)

@pytest.fixture
def complex_dummy_trajectory():
    cell_ids = [f"cell{i}" for i in range(5)]
    expr = pd.DataFrame(
        np.random.rand(5, 5),
        index=cell_ids,
        columns=[f"gene{i}" for i in range(5)]
    )
    milestone_percentages = pd.DataFrame({
        "cell_id": cell_ids * 3,
        "milestone_id": ["m1"]*5 + ["m2"]*5 + ["m3"]*5,
        "percentage": np.random.rand(15)
    })
    return {
        "cell_ids": cell_ids,
        "expression": expr,
        "milestone_percentages": milestone_percentages,
        "milestone_ids": ["m1", "m2", "m3"],
        "pydynwrap:with_trajectory": True
    }

def test_complex_calculate_overall_feature_importance(monkeypatch, complex_dummy_trajectory):
    # 替换为新的模拟函数（显式补全所有组合）
    monkeypatch.setattr(
        "pydynverse.feature.calculate_overall_feature_importance.calculate_milestone_feature_importance",
        dummy_milestone_feature_importance_complex
    )
    
    # 调用被测函数
    dummy_fi_method = {"fun": lambda X, y, verbose: {}}
    result = calculate_overall_feature_importance(
        trajectory=complex_dummy_trajectory,
        expression_source="expression",
        fi_method=dummy_fi_method,
        verbose=False
    )
    
    # 预期结果：每个特征在所有3个里程碑中的平均重要性（包括显式补零）
    expected_data = {
        "feature_id": ["gene0", "gene1", "gene2", "gene3", "gene4"],
        "importance": [
            (0.9 + 0.7 + 0.5) / 3,  # gene0
            (0.6 + 0.8 + 0) / 3,     # gene1
            (0.0 + 0 + 0.3) / 3,     # gene2
            (0 + 0.4 + 0) / 3,       # gene3
            (0 + 0 + 0.6) / 3        # gene4
        ]
    }
    expected = pd.DataFrame(expected_data)
    
    # 排序并断言
    result_sorted = result.sort_values("feature_id").reset_index(drop=True)
    expected_sorted = expected.sort_values("feature_id").reset_index(drop=True)
    pd.testing.assert_frame_equal(result_sorted, expected_sorted, atol=1e-9)

if __name__ == "__main__":
    pytest.main(["-v", __file__])