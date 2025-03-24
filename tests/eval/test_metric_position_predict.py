import numpy as np
import pandas as pd
import pytest
from sklearn.linear_model import LinearRegression  
from sklearn.ensemble import RandomForestRegressor  
from pydynverse.util.expand_matrix import expand_matrix
from pydynverse.eval.metric_position_predict import calculate_position_predict  

# 构造一个简单的 dummy 数据集和 prediction 用于测试
@pytest.fixture
def dummy_dataset():
    # 模拟一个 dataset，包含 cell_ids 与 milestone_percentages
    cell_ids = ["c1", "c2", "c3", "c4"]
    # milestone_percentages 为长格式 DataFrame，每个细胞在多个里程碑上的百分比
    milestone_percentages = pd.DataFrame({
        "cell_id": ["c1", "c1", "c2", "c2", "c3", "c3", "c4", "c4"],
        "milestone_id": ["m1", "m2", "m1", "m2", "m1", "m2", "m1", "m2"],
        "percentage": [0.8, 0.2, 0.5, 0.5, 1.0, 0.0, 0.3, 0.7]
    })
    return {
        "cell_ids": cell_ids,
        "milestone_percentages": milestone_percentages
    }

@pytest.fixture
def dummy_prediction():
    # 构造一个 prediction，其 milestone_percentages 数据与 dataset 略有不同
    cell_ids = ["c1", "c2", "c3", "c4"]
    milestone_percentages = pd.DataFrame({
        "cell_id": ["c1", "c1", "c2", "c2", "c3", "c3", "c4", "c4"],
        "milestone_id": ["m1", "m2", "m1", "m2", "m1", "m2", "m1", "m2"],
        # 模拟预测值与真实值存在一定差异，例如将 dataset 的值乘以一个因子，再加上噪声
        "percentage": [0.85, 0.15, 0.55, 0.45, 0.95, 0.05, 0.35, 0.65]
    })
    return {
        "milestone_percentages": milestone_percentages
    }

@pytest.fixture
def dummy_dataset_with_ids(dummy_dataset):
    # 为了调用 calculate_position_predict，dataset 需要包含 "cell_ids" 键和表达矩阵等
    # 这里我们基于 dummy_dataset 补充一个简单的表达矩阵，
    # 表达矩阵的行索引必须与 cell_ids 对应，这里我们构造一个 4 x 3 的矩阵
    cell_ids = dummy_dataset["cell_ids"]
    expression = pd.DataFrame({
        "gene1": [1, 2, 3, 4],
        "gene2": [2, 3, 4, 5],
        "gene3": [3, 4, 5, 6],
    }, index=cell_ids)
    # 将表达矩阵存入 dataset 字典中
    dataset = dummy_dataset.copy()
    dataset["expression"] = expression
    return dataset

@pytest.fixture
def dummy_prediction_with_ids(dummy_prediction, dummy_dataset):
    # 补充 prediction 中的 cell_ids 和表达矩阵（这里直接用 dataset 的表达矩阵进行测试）
    cell_ids = dummy_dataset["cell_ids"]
    # 可采用与 dataset 略有差异的表达矩阵，模拟预测的变化，比如将数值稍作放大
    expression = pd.DataFrame({
        "gene1": [1.1, 2.1, 3.1, 4.1],
        "gene2": [2.1, 3.1, 4.1, 5.1],
        "gene3": [3.1, 4.1, 5.1, 6.1],
    }, index=cell_ids)
    prediction = dummy_prediction.copy()
    prediction["cell_ids"] = cell_ids
    prediction["expression"] = expression
    return prediction

def test_calculate_position_predict_normal(dummy_dataset_with_ids, dummy_prediction_with_ids):
    # 测试正常情况下的指标计算
    output = calculate_position_predict(
        dataset=dummy_dataset_with_ids,
        prediction=dummy_prediction_with_ids,
        metrics=["rf_mse", "rf_rsq", "lm_mse", "lm_rsq"]
    )
    # 检查输出为字典，并且 summary 中至少包含 rf_mse, rf_rsq, lm_mse, lm_rsq
    assert isinstance(output, dict)
    summary = output.get("summary", {})
    for key in ["rf_mse", "rf_rsq", "lm_mse", "lm_rsq"]:
        assert key in summary, f"Missing key {key} in summary"
        # 指标值应为数值且非负
        assert isinstance(summary[key], (int, float))
        assert summary[key] >= 0

def test_calculate_position_predict_invalid(dummy_dataset_with_ids):
    # 当 prediction 为 None 或 cell_id 数量不足时，应返回默认 summary
    output = calculate_position_predict(
        dataset=dummy_dataset_with_ids,
        prediction=None,
        metrics=["rf_mse", "rf_rsq", "lm_mse", "lm_rsq"]
    )
    summary = output.get("summary", {})
    # 此时应返回 rf_mse 与 lm_mse 为 baseline（计算自 gold_milenet_m），其他为 0
    assert "rf_mse" in summary and "lm_mse" in summary
    # 这里简单检查非空即可
    assert summary["rf_mse"] is not None
    assert summary["lm_mse"] is not None

if __name__ == "__main__":
    pytest.main(["-v", __file__])
