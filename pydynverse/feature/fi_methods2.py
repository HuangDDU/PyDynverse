import inspect
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
# TODO:新版本（version0.0.2：2版本文件命令会在后边加一个2，依赖接口类似） 

def apply_function_params(params: dict, nrow, ncol):
    """
    遍历字典中的每个参数，如果参数是可调用对象且其签名正好要求 nrow 和 ncol，
    则用提供的 nrow, ncol 调用该函数，并将返回值替换原参数值。
    """
    new_params = {}
    for key, value in params.items():
        if callable(value):
            sig = inspect.signature(value)
            param_names = list(sig.parameters.keys())
            if len(param_names) == 2 and set(param_names) == {"nrow", "ncol"}:
                new_params[key] = value(nrow=nrow, ncol=ncol)
            else:
                new_params[key] = value
        else:
            new_params[key] = value
    return new_params

def fi_ranger_rf(num_trees, mtry, sample_fraction, min_node_size, **kwargs):
    """
    构造一个基于 scikit-learn RandomForestClassifier 的特征重要性函数。
    """
    default_params = {
        "n_jobs": 1,
        "min_samples_leaf": min_node_size,
        # 模拟 R 代码中的默认设置
        "importance": "impurity", 
        "write_forest": False
    }
    params = {**default_params, **kwargs}
    
    def fi_function(X, y, verbose=False):
        # 如果 X 不是 DataFrame，则转换
        if not isinstance(X, pd.DataFrame):
            X = pd.DataFrame(X)
        nrow, ncol = X.shape
        
        # 计算 mtry 和 sample_fraction
        max_features = mtry(nrow=nrow, ncol=ncol)
        fraction = sample_fraction(nrow=nrow, ncol=ncol)
        max_samples = fraction if fraction < 1 else None
        
        # 构造数据：构造一个新的 DataFrame，其中第一列为目标变量 "PREDICT"，后续为特征 X
        data = X.copy()
        data.insert(0, "PREDICT", y)
        
        rf = RandomForestClassifier(
            n_estimators=num_trees,
            max_features=max_features,
            min_samples_leaf=params.get("min_samples_leaf"),
            n_jobs=params.get("n_jobs"),
            max_samples=max_samples,
            bootstrap=True if max_samples is not None else False,
            random_state=42
        )
        if verbose:
            print("训练随机森林，参数：", rf.get_params())
        # 使用除 PREDICT 之外的所有列作为特征进行训练
        rf.fit(data.drop("PREDICT", axis=1), data["PREDICT"])
        importance = rf.feature_importances_
        return dict(zip(X.columns, importance))
    
    return {"fun": fi_function}

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
    模拟 R 中 caret 接口的特征重要性函数（仅支持 'rf'）。
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
    
    return {"fun": fi_function}

def fi_ranger_rf_tiny(num_trees=100, num_variables_per_split=50, num_samples_per_tree=250, min_node_size=20, **kwargs):
    """
    fi_ranger_rf 的轻量级小型版本，默认树的数量较少。
    """
    def mtry(nrow, ncol):
        return min(num_variables_per_split, ncol)
    def sample_fraction(nrow, ncol):
        return min(num_samples_per_tree / nrow, 1)
    return fi_ranger_rf(num_trees, mtry, sample_fraction, min_node_size, **kwargs)
