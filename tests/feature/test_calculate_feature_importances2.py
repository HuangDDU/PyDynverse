import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csc_matrix
from pydynverse.feature.calculate_feature_importances2 import calculate_feature_importances
from pydynverse.feature.fi_methods2 import fi_ranger_rf_lite

# 构造一个 dummy fi_method 用于测试对比（返回每列样本方差）
def dummy_fi_method():
    def fun(X, y, verbose=False):
        # 返回 X 各列样本方差
        variances = X.var(ddof=1)
        return variances.to_dict()
    return {"fun": fun}

@pytest.fixture
def X_df():
    # 构造 100x5 的 DataFrame
    np.random.seed(42)
    data = np.random.randn(100, 5)
    return pd.DataFrame(data, columns=[f"feature_{i}" for i in range(5)])

@pytest.fixture
def Y_df():
    # 构造 100x3 的 DataFrame，其中所有列均为离散类别（例如二分类）
    np.random.seed(0)
    col1 = np.random.randint(0, 2, 100)
    col2 = np.random.randint(0, 2, 100)  # 修改为离散
    col3 = np.random.randint(0, 2, 100)  # 修改为离散
    return pd.DataFrame({"pred1": col1, "pred2": col2, "pred3": col3})

def test_calculate_feature_importances_normal(X_df, Y_df):
    # 使用默认的 fi_method（fi_ranger_rf_lite()）测试
    result = calculate_feature_importances(X_df, Y_df, fi_method=fi_ranger_rf_lite(), verbose=False)
    # 检查返回 DataFrame 是否包含 predictor_id, feature_id, importance 三列
    assert isinstance(result, pd.DataFrame)
    expected_cols = {"predictor_id", "feature_id", "importance"}
    assert set(result.columns) == expected_cols, f"预期列 {expected_cols}, 得到 {set(result.columns)}"
    # 对于 Y 中常数列情况测试：这里 Y_df 不存在全常数列
    # 检查结果是否按 importance 降序排列
    imp_vals = result["importance"].values
    assert np.all(np.diff(imp_vals) <= 0), "结果未按 importance 降序排列"

def test_calculate_feature_importances_Y_array(X_df):
    # 当 Y 为 numpy 数组时，转换后其列名为整数 0,1,2
    Y = np.random.randint(0, 2, size=(100, 3))
    result = calculate_feature_importances(X_df, Y, fi_method=dummy_fi_method(), verbose=False)
    expected_predictors = {0, 1, 2}  # 这里预期为整数
    assert set(result["predictor_id"].unique()) == expected_predictors

def test_calculate_feature_importances_Y_sparse(X_df):
    # 当 Y 为稀疏矩阵时
    Y_dense = np.random.randint(0, 2, size=(100, 3))
    Y_sparse = csc_matrix(Y_dense)
    result = calculate_feature_importances(X_df, Y_sparse, fi_method=dummy_fi_method(), verbose=False)
    expected_predictors = {0, 1, 2}  # 预期为整数
    assert set(result["predictor_id"].unique()) == expected_predictors

if __name__ == "__main__":
    pytest.main(["-v", __file__])
