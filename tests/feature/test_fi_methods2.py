import numpy as np
import pandas as pd
import pytest
from sklearn.datasets import make_classification

from pydynverse.feature.fi_methods2 import (
    apply_function_params,
    fi_ranger_rf_lite,
    fi_caret,
    fi_ranger_rf_tiny,
)

# 构造一个简单的分类数据集，用于测试特征重要性计算
@pytest.fixture
def classification_data():
    # 生成一个二分类数据集，100个样本，10个特征
    X, y = make_classification(n_samples=100, n_features=10, n_informative=5, random_state=42)
    # 转换为 DataFrame，并设置列名
    X_df = pd.DataFrame(X, columns=[f"feature_{i}" for i in range(10)])
    # 将 y 转换为离散标签（分类问题）
    y_series = pd.Series(y, name="target")
    return X_df, y_series

def test_apply_function_params():
    # 定义一个测试函数，返回 nrow * ncol
    def func(nrow, ncol):
        return nrow * ncol
    params = {"a": 1, "b": func, "c": lambda x: x}
    result = apply_function_params(params, nrow=4, ncol=5)
    # 参数 a 应不变，参数 b 应为20，参数 c 不匹配签名保持原函数
    assert result["a"] == 1
    assert result["b"] == 20
    assert callable(result["c"])

def test_fi_ranger_rf_lite(classification_data):
    X, y = classification_data
    # 获取特征重要性函数
    fi_obj = fi_ranger_rf_lite(num_trees=100, num_variables_per_split=5, num_samples_per_tree=50, min_node_size=2)
    assert "fun" in fi_obj
    fi_fun = fi_obj["fun"]
    importance = fi_fun(X, y, verbose=True)
    # 检查返回结果是字典，键为 X 的列名，且所有值为非负数
    assert isinstance(importance, dict)
    assert set(importance.keys()) == set(X.columns)
    for imp in importance.values():
        assert imp >= 0

def test_fi_caret(classification_data):
    X, y = classification_data
    # 测试 fi_caret 仅支持 "rf"
    fi_obj = fi_caret("rf")
    fi_fun = fi_obj["fun"]
    importance = fi_fun(X, y, verbose=True)
    # 检查返回结果结构与 fi_ranger_rf_lite 一致
    assert isinstance(importance, dict)
    assert set(importance.keys()) == set(X.columns)

def test_fi_ranger_rf_tiny(classification_data):
    X, y = classification_data
    fi_obj = fi_ranger_rf_tiny()
    fi_fun = fi_obj["fun"]
    importance = fi_fun(X, y, verbose=True)
    # 检查返回结果
    assert isinstance(importance, dict)
    assert set(importance.keys()) == set(X.columns)

def test_integration_all_methods(classification_data):
    X, y = classification_data
    # 分别调用三种方法
    methods = {
        "lite": fi_ranger_rf_lite(),
        "caret": fi_caret("rf"),
        "tiny": fi_ranger_rf_tiny()
    }
    results = {}
    for key, fi_obj in methods.items():
        fi_fun = fi_obj["fun"]
        imp = fi_fun(X, y)
        results[key] = imp
    # 检查三个方法的返回字典都包含相同的键（即所有特征名）
    expected_keys = set(X.columns)
    for res in results.values():
        assert set(res.keys()) == expected_keys
    # 如果数据集足够稳定，三个方法返回的特征重要性数值应大致接近（差异可能由随机森林差异引起，但不宜太大）
    for feat in expected_keys:
        vals = [results[m][feat] for m in results]
        mean_val = np.mean(vals)
        std_val = np.std(vals)
        if mean_val < 1e-4:
            # 对于不重要的特征，所有方法的结果都应很低
            assert max(vals) < 0.05, f"Feature {feat} 不重要，所有方法的值应接近0，但最大值为 {max(vals)}"
        else:
            cv = std_val / mean_val  # 计算变异系数
            assert cv < 2.0, f"Feature {feat} 方法间变异系数过高: {cv}"

if __name__ == "__main__":
    pytest.main(["-v", __file__])
