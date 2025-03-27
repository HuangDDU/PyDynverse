import pydynverse
from pydynverse.feature.fi_methods import *
import pytest

def test_apply_function_params():
    # 定义一个符合要求的函数，返回 nrow 与 ncol 的乘积
    def func(nrow, ncol):
        return nrow * ncol
    # 定义一个不符合签名的 lambda 函数
    params = {
        'a': 1,
        'b': func,
        'c': lambda x: x  # 参数不匹配时应保持不变
    }
    result = apply_function_params(params, nrow=4, ncol=5)
    assert result['a'] == 1, "非函数参数应保持不变"
    assert result['b'] == 20, "函数参数应被调用，4*5=20"
    # 检查 'c' 仍为函数（未被调用）
    assert callable(result['c']), "签名不匹配的 lambda 应保持不变"

def test_fi_ranger_rf_lite():
    np.random.seed(42)
    # 生成 100x10 的测试数据，列名为 feature0, feature1, ..., feature9
    X = pd.DataFrame(np.random.randn(100, 10), columns=[f'feature{i}' for i in range(10)])
    # 生成样本的分类标签
    y = np.random.randint(0, 2, size=100)
    model_wrapper = fi_ranger_rf_lite()
    assert 'fun' in model_wrapper, "返回的包装字典中应包含 'fun' 键"
    fi_fun = model_wrapper['fun']
    importance = fi_fun(X, y, verbose=True)
    print(f'importance is {importance}')


def test_fi_ranger_rf_tiny():
    np.random.seed(42)#固定种子
    # 生成 50x5 的测试数据
    X = pd.DataFrame(np.random.randn(50, 5), columns=[f'feat{i}' for i in range(5)])
    y = np.random.randint(0, 2, size=50)
    model_wrapper = fi_ranger_rf_tiny()
    fi_fun = model_wrapper['fun']
    importance = fi_fun(X, y)
    print(f'importance is {importance}')

def test_fi_caret():
    np.random.seed(42)#固定种子
    # 生成 80x8 的测试数据
    X = pd.DataFrame(np.random.randn(80, 8), columns=[f'col{i}' for i in range(8)])
    y = np.random.randint(0, 2, size=80)
    model_wrapper = fi_caret("rf")
    fi_fun = model_wrapper['fun']
    importance = fi_fun(X, y, verbose=True)
    print(f'importance is {importance}')


if __name__ == "__main__":
    pytest.main(["-v", __file__])
