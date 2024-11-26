import numpy as np
import scipy.sparse as sp

from .._logging import logger
from .wrap_data import wrap_data


def dgC2csc(dgC_obj):
    # 提取列压缩矩阵的行、列、值向量
    row_indices = np.array(dgC_obj.slots["i"])
    col_indices = np.array(dgC_obj.slots["p"])
    value_array = np.array(dgC_obj.slots["x"])
    shape = np.array(dgC_obj.slots["Dim"])

    csc = sp.csc_matrix((value_array, row_indices, col_indices), shape=shape)

    return csc


def add_expression(data, counts, expression):

    counts_csc = dgC2csc(counts)
    expression_csc = dgC2csc(expression)

    data["counts"] = counts_csc
    data["expression"] = expression_csc

    # 同时修改内部的adata
    adata = data["adata"]
    adata.layers["counts"] = counts_csc
    adata.X = counts_csc

    return data


def is_wrapper_with_expression(dataset):
    return True


def get_expression(dataset, expression_source="expression"):
    # 提取表达矩阵
    expression_source_type = type(expression_source)
    if expression_source_type == str:
        # 根据键提取表达矩阵值
        expression = dataset[expression_source]
    elif expression_source_type == np.ndarray or expression_source_type == sp._csc.csc_matrix:
        # 直接作为表达矩阵
        expression = expression_source
    return expression


def wrap_expression(expression,
                    counts,
                    id=None,
                    cell_info=None,
                    feature_info=None,
                    expression_future=None):
    # 封装让rpy2的expression, counts稀疏矩阵为dataset字典的一项

    # 解析R对象，获得obs属性
    cell_ids, feature_ids = counts.slots["Dimnames"]
    cell_ids = np.array(cell_ids)
    feature_ids = np.array(feature_ids)
    logger.debug("Resolving cell_ids and feature_ids")

    # 创建AnnData对象, 初始表达0
    dataset = wrap_data(
        cell_ids,
        feature_ids,
        id=id,
        cell_info=cell_info,
        feature_info=feature_info,
    )
    logger.debug(f"Dataset created: {dataset}")

    # 添加表达矩阵到AnnData对象
    dataset = add_expression(dataset, counts, expression)
    logger.debug(f"Andata counts and expression added")

    return dataset
