import numpy as np
import pandas as pd
import pytest
from pydynverse.feature.calculate_milestone_feature_importance import calculate_milestone_feature_importance

# 定义一个简单的 dummy 特征重要性方法：
def dummy_fi_method():
    """
    对于给定的 X, y，返回各列方差作为重要性分数。
    """
    def fun(X, y, verbose=False):
        # 使用每列的方差作为“重要性”，返回字典格式
        return X.var().to_dict()
    return {'fun': fun}

# 为保证测试时 get_expression 与 is_wrapper_with_trajectory 能够正常调用，
# 定义 dummy 实现，并将其放入全局名称空间，供 calculate_milestone_feature_importance 使用。
def dummy_get_expression(trajectory, expression_source="expression"):
    # 假设 trajectory[expression_source] 即为表达矩阵
    return trajectory[expression_source]

def dummy_is_wrapper_with_trajectory(trajectory):
    return trajectory.get("pydynwrap:with_trajectory", False)

# 将全局名称绑定（如果 calculate_milestone_feature_importance 中直接使用了这两个名称）
globals()["get_expression"] = dummy_get_expression
globals()["is_wrapper_with_trajectory"] = dummy_is_wrapper_with_trajectory

def test_calculate_milestone_feature_importance():
    # 构造模拟的轨迹数据
    cell_ids = ["cell1", "cell2", "cell3", "cell4"]
    # 构造表达矩阵：4 个细胞，5 个基因
    expr = pd.DataFrame(np.random.rand(4, 5), index=cell_ids, columns=[f"gene{i}" for i in range(5)])
    
    # 构造 milestone_percentages 的长格式 DataFrame，
    # 每个细胞有两个里程碑 "m1" 和 "m2"，给出各自的百分比
    data = [
        {"cell_id": "cell1", "milestone_id": "m1", "percentage": 0.8},
        {"cell_id": "cell1", "milestone_id": "m2", "percentage": 0.2},
        {"cell_id": "cell2", "milestone_id": "m1", "percentage": 0.5},
        {"cell_id": "cell2", "milestone_id": "m2", "percentage": 0.5},
        {"cell_id": "cell3", "milestone_id": "m1", "percentage": 1.0},
        {"cell_id": "cell3", "milestone_id": "m2", "percentage": 0.0},
        {"cell_id": "cell4", "milestone_id": "m1", "percentage": 0.3},
        {"cell_id": "cell4", "milestone_id": "m2", "percentage": 0.7},
    ]
    milestone_percentages = pd.DataFrame(data)
    
    # 构造 trajectory 字典
    trajectory = {
        "cell_ids": cell_ids,
        "milestone_percentages": milestone_percentages,
        "milestone_ids": ["m1", "m2"],
        "pydynwrap:with_trajectory": True,
        "expression": expr
    }
    
    # 调用 calculate_milestone_feature_importance，使用 dummy_fi_method 计算特征重要性
    result = calculate_milestone_feature_importance(
        trajectory, 
        expression_source="expression", 
        milestones_oi=None,  # 使用所有里程碑
        fi_method=dummy_fi_method(), 
        verbose=False
    )
    
    # 检查返回结果：应包含 "milestone_id", "feature_id", "importance" 三列
    assert isinstance(result, pd.DataFrame)
    expected_cols = {"milestone_id", "feature_id", "importance"}
    assert set(result.columns) == expected_cols, f"Expected columns {expected_cols}, got {set(result.columns)}"
    
    # 检查里程碑标签只包含 "m1" 和 "m2"
    milestone_values = set(result["milestone_id"].unique())
    assert milestone_values.issubset({"m1", "m2"})
    
    # 结果不应为空
    assert not result.empty

if __name__ == "__main__":
    pytest.main(["-v", __file__])
