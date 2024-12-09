from typing import Optional, Literal, get_args
import inspect

import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import pairwise_distances


def calculate_distance(
    x: pd.DataFrame,
    y: pd.DataFrame = None,
    method: Optional[Literal["manhattan", "euclidean", "cosine", "correlation"]] = None,
    margin: int = 1,  # 1标识行为样本，2标识列为样本
    diag: bool = False,
    drop0: bool = False
):
    input_ = _process_input_matrices(x, y, margin)
    x = input_["x"]
    y = input_["y"]

    # 这里与R逻辑不同，不用手动转换相似性与距离，sklearn内部会自动进行的转换
    dis = pairwise_distances(x, y, metric=method)  # 获得的是array结果

    if y is None:
        np.fill_diagonal(dis, 0)

    dis = pd.DataFrame(dis, index=x.index, columns=y.index)  # 转化为DataFrame

    return dis


def list_similarity_methods():
    signature = inspect.signature(calculate_distance)  # 获取函数签名
    method_annotation = signature.parameters["method"].annotation
    return list(get_args(get_args(method_annotation)[0]))


def _process_input_matrices(x, y, margin=1):

    # 维度转化为 样本数*特征数
    if margin == 2:
        x = x.T
        y = y.T

    input_ = {"x": x, "y": y}

    return input_
