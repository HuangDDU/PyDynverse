import numpy as np
import pandas as pd
from scipy.sparse import issparse
from pydynverse.feature.fi_methods2 import fi_ranger_rf_lite
#TODO：version0.0.2

def calculate_feature_importances(X, Y, fi_method=fi_ranger_rf_lite(), verbose=False):
    """
    计算特征重要性，返回一个 DataFrame，包含 predictor_id, feature_id 和 importance 三列。
    
    Parameters:
      X : pd.DataFrame or array-like
          包含特征的矩阵，特征为列。
      Y : pd.DataFrame or array-like
          包含预测变量的矩阵，其行数应与 X 相同。如果不是 DataFrame，
          则如果为稀疏矩阵则转换为稠密矩阵，再转换为 DataFrame。
      fi_method : dict
          特征重要性方法，字典中键 "fun" 对应一个函数，签名为 fun(X, y, verbose)。
          默认使用 fi_ranger_rf_lite() 的返回值。
      verbose : bool, optional
          是否输出详细信息。
    
    Returns:
      pd.DataFrame : 包含 predictor_id, feature_id, importance 三列。
          predictor_id 为 Y 的列名（如果 Y 原来是 DataFrame，则使用其列名，
          如果是数组，则使用整数），feature_id 为 X 的列名（对应 fi_method 计算结果中的键），
          importance 为计算得到的重要性数值。最终结果按 importance 降序排列。
    """
    # 如果 Y 不是 DataFrame，则转换
    if not isinstance(Y, pd.DataFrame):
        if issparse(Y):
            Y = pd.DataFrame(Y.toarray())
        else:
            Y = pd.DataFrame(Y)
    
    # 如果 X 是稀疏矩阵，则转换为稠密矩阵
    if issparse(X):
        X = pd.DataFrame(X.toarray())
    elif not isinstance(X, pd.DataFrame):
        X = pd.DataFrame(X)
    
    result_list = []
    # 遍历 Y 的每一列（每个预测变量）
    for i, predictor in enumerate(Y.columns):
        if verbose:
            print(f"Generating forest {i+1}/{len(Y.columns)}")
        
        y = Y[predictor]
        # 如果 y 类型为字符或布尔型，则转换为类别变量
        if y.dtype == object or y.dtype == bool:
            y = pd.Categorical(y)
        
        # 如果 y 的所有值相同，则 importance 全部为 0
        if len(y.unique()) == 1:
            importance_dict = {feat: 0 for feat in X.columns}
        else:
            importance_dict = fi_method["fun"](X, y, verbose=verbose)
        
        # 构造当前 predictor 的 DataFrame
        df = pd.DataFrame({
            "predictor_id": predictor,
            "feature_id": list(importance_dict.keys()),
            "importance": list(importance_dict.values())
        })
        result_list.append(df)
    
    result_df = pd.concat(result_list, ignore_index=True)
    result_df = result_df.sort_values(by="importance", ascending=False).reset_index(drop=True)
    return result_df
