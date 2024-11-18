import numpy as np

from .._logging import logger
from .wrap_data import wrap_data


def dgC2csc(dgC_obj):
    from scipy.sparse import csc_matrix

    row_indices = np.array(dgC_obj.slots["i"])
    col_indices = np.array(dgC_obj.slots["p"])
    value_array = np.array(dgC_obj.slots["x"])
    shape = np.array(dgC_obj.slots["Dim"])

    csc = csc_matrix((value_array, row_indices, col_indices), shape=shape)

    return csc


def add_expression(adata, counts, expression):
    adata.layers["counts"] = dgC2csc(counts)
    adata.X = dgC2csc(expression)

    return adata


def wrap_expression(expression, counts, id=None, cell_info=None, feature_info=None, expression_future=None):
    # 解析R对象，获得obs属性
    cell_ids, feature_ids = counts.slots["Dimnames"]
    cell_ids = np.array(cell_ids)
    feature_ids = np.array(feature_ids)
    logger.debug("Resolving cell_ids and feature_ids")

    # 创建AnnData对象, 初始表达0
    adata = wrap_data(
        cell_ids,
        feature_ids,
        id=id,
        cell_info=cell_info,
        feature_info=feature_info,
    )
    logger.debug(f"AnnData created: {adata}")
    

    # 添加表达矩阵到AnnData对象
    adata = add_expression(adata, counts, expression)
    logger.debug(f"Andata counts and expression added")

    return adata
