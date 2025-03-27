import inspect
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
#可能存在与R逻辑不一致问题
#TODO:(初步实现)version0.0.1

def apply_function_params(params:dict, nrow, ncol):
    """
    遍历字典中的每个参数，如果参数是可调用对象且其签名正好要求 nrow 和 ncol，
    则用提供的 nrow, ncol 调用该函数，并将返回值替换原参数值。
    """
    new_params = {}
    for key, value in params.items():
        if callable(value):
            # 获取函数签名
            sig = inspect.signature(value)
            param_names = list(sig.parameters.keys())
            # 检查是否正好包含 nrow 和 ncol
            if len(param_names) == 2 and set(param_names) == {'nrow', 'ncol'}:
                new_params[key] = value(nrow=nrow, ncol=ncol)
            else:
                new_params[key] = value
        else:
            new_params[key] = value
    return new_params


def fi_ranger_rf(num_trees, mtry, sample_fraction, min_node_size, **kwargs):
    #这里的mtry是否需要指定默认？
    """
     构造一个基于 scikit-learn RandomForestClassifier 的特征重要性函数
    Args:
        num_trees (_type_): 随机森林的树数量
        mtry (_type_): 一个函数，根据 (nrow, ncol) 返回每次分裂时考虑的特征数
        sample_fraction (_type_):  一个函数，根据 (nrow, ncol) 返回样本抽样比例
        min_node_size (_type_):  叶子节点最小样本数
        **kwards:其他额外参数
    Returns:
        _type_: _description_
    """
    # 默认参数设置，类似 R 中 default_params
    default_params = {
        'n_jobs': 1,  # 对应 num.threads=1
        'min_samples_leaf': min_node_size
        # 其他参数可以在 kwargs 中传入
    }
    # 合并额外参数，优先 kwargs（类似 R 中 list_modify）
    params = {**default_params, **kwargs}

    def fi_function(X, y, verbose=False):
        # 保证 X 为二维数组或 DataFrame
        if isinstance(X, pd.DataFrame):
            nrow, ncol = X.shape
        else:
            X = np.asarray(X)
            nrow, ncol = X.shape

        # 计算 mtry 和 sample_fraction 对应的数值
        max_features = mtry(nrow=nrow, ncol=ncol)
        fraction = sample_fraction(nrow=nrow, ncol=ncol)
        # scikit-learn 的 RandomForestClassifier 提供 max_samples（0-1 时表示比例）
        max_samples = fraction if fraction < 1 else None

        # 创建随机森林模型
        # scikit-learn 中使用 max_features 参数来控制每个拆分时考虑的特征数
        rf = RandomForestClassifier(
            n_estimators=num_trees,
            max_features=max_features,
            min_samples_leaf=params.get('min_samples_leaf'),
            n_jobs=params.get('n_jobs'),
            max_samples=max_samples,
            bootstrap=True if max_samples is not None else False,
            random_state=42  # 固定随机种子以便复现
        )
        if verbose:
            print("训练随机森林，参数：", rf.get_params())
        # 拟合模型
        rf.fit(X, y)
        # 获取特征重要性
        importance = rf.feature_importances_
        # 如果 X 为 DataFrame，则返回 {特征名: 重要性} 的字典，否则返回 {索引: 重要性}
        if isinstance(X, pd.DataFrame):
            return dict(zip(X.columns, importance))
        else:
            return dict(enumerate(importance))
    # 返回一个字典，其中 'fun' 键对应特征重要性计算的函数（类似 R 中返回的 list(fun = function(...){ ... })）
    return {'fun': fi_function}


def fi_ranger_rf_lite(num_trees=2000, num_variables_per_split=50, num_samples_per_tree=250, min_node_size=20, **kwargs):
    """
    轻量级版本，封装了 fi_ranger_rf 的参数，并提供默认值。
    
    其中：
     mtry 函数：返回 min(num_variables_per_split, ncol)
     sample_fraction 函数：返回 min(num_samples_per_tree / nrow, 1)
    """
    def mtry(nrow, ncol):
        return min(num_variables_per_split, ncol)
    def sample_fraction(nrow, ncol):
        return min(num_samples_per_tree / nrow, 1)
    return fi_ranger_rf(num_trees, mtry, sample_fraction, min_node_size, **kwargs)


def fi_caret(caret_method, **kwargs):
    """
    模拟 R 中 caret 接口的特征重要性函数。
    由于 Python 中没有 caret 包，这里简单地将 'rf' 映射到随机森林。
    
    参数：
     caret_method: 目前仅支持 "rf"，否则抛出异常
     kwargs: 其他额外参数传递给 RandomForestClassifier
    """
    if caret_method != "rf":
        raise ValueError("Invalid method. Only 'rf' is supported in this demo.")
    
    def fi_function(X, y, verbose=False):
        rf = RandomForestClassifier(random_state=42, **kwargs)
        rf.fit(X, y)
        importance = rf.feature_importances_
        if isinstance(X, pd.DataFrame):
            return dict(zip(X.columns, importance))
        else:
            return dict(enumerate(importance))
    return {'fun': fi_function}


def fi_ranger_rf_tiny(num_trees=100, num_variables_per_split=50, num_samples_per_tree=250, min_node_size=20, **kwargs):
    """
    fi_ranger_rf 的轻量级小型版本，默认树的数量较少。
    参数设定同 fi_ranger_rf_lite，只是默认 num_trees 由 2000 改为 100。
    """
    def mtry(nrow, ncol):
        return min(num_variables_per_split, ncol)
    def sample_fraction(nrow, ncol):
        return min(num_samples_per_tree / nrow, 1)
    return fi_ranger_rf(num_trees, mtry, sample_fraction, min_node_size, **kwargs)
