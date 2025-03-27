from pydynverse.feature.calculate_feature_importances import calculate_feature_importances
import pytest
import numpy as np
import pandas as pd
from tqdm import tqdm

def dummy_fi_method():
    """
    这里定义一个简单的 dummy 方法：
    对于给定的 X, y，返回每个特征的方差作为重要性分数。
    这样既简单又能体现不同特征的重要性差异。
    """
    def fun(X, y, verbose=False):
        importance = X.var().to_dict()  # 直接使用各列的方差
        return importance
    return {'fun': fun}

# 测试 calculate_feature_importances 的函数

def test_calculate_feature_importances():

    np.random.seed(42)
    # 生成一个 25 x 10 的随机数据集 X，特征名为 feature0, feature1, ..., feature9
    X = pd.DataFrame(np.random.rand(25, 10), columns=[f'feature{i}' for i in range(10)])
    # 生成 Y，只有一列，所有值均为 1（常量预测器）
    Y_const = pd.DataFrame(np.ones((25, 1)), columns=['const'])
    # 使用 dummy 方法计算特征重要性
    result_const = calculate_feature_importances(X, Y_const, fi_method=dummy_fi_method(), verbose=True)
    # 对于常量预测器，重要性应该全部为 0
    assert (result_const['importance'] == 0).all(), "常量预测器的所有重要性应为 0"
    
    # 生成 Y，只有一列，值为随机数（非常量预测器）
    Y_nonconst = pd.DataFrame(np.random.rand(25, 1), columns=['non_const'])
    result_nonconst = calculate_feature_importances(X, Y_nonconst, fi_method=dummy_fi_method(), verbose=True)
    # 对于非常量预测器，dummy 方法返回的各特征方差应该大于 0（总体上非零）
    assert result_nonconst['importance'].sum() > 0, "非常量预测器应返回非零重要性值"
    
    # 构造 Y 包含两个预测器，一列为随机数（非常量），另一列为常量
    Y_multi = pd.DataFrame({
        'pred1': np.random.rand(25),
        'pred2': np.ones(25)  # 常量预测器
    })
    result_multi = calculate_feature_importances(X, Y_multi, fi_method=dummy_fi_method(), verbose=True)
    # 对于 pred1，应返回非零重要性；而对于 pred2，应返回全零
    pred1_imp = result_multi[result_multi['predictor_id'] == 'pred1']['importance']
    pred2_imp = result_multi[result_multi['predictor_id'] == 'pred2']['importance']
    assert pred1_imp.sum() > 0, "非常量预测器 pred1 应返回非零重要性值"
    assert (pred2_imp == 0).all(), "常量预测器 pred2 应返回所有重要性为 0"



if __name__ == "__main__":
    pytest.main(["-v", __file__])

