import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from .._logging import logger
from .wrap_data import wrap_data


def add_expression(data, counts, expression):

    # 传入的就是csc矩阵
    data["counts"] = counts
    data["expression"] = expression

    feature_ids = data["feature_ids"] 
    if feature_ids is None:
        # 之前没有提供feature_ids, 此处根据序号命名
        feature_ids = [f"feature_{i}" for i in range(data["counts"].shape[1])]
        feature_info = pd.DataFrame(index=feature_ids)
        data["feature_ids"] = feature_ids
        data["feature_info"] = feature_info

    # 创建adata
    adata = ad.AnnData(
        X=counts,
        obs=data["cell_info"],
        var=data["feature_info"],
    )
    adata.layers["counts"] = counts
    adata.layers["expression"] = expression

    data["adata"] = adata

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
    # 从expression矩阵的行名和列名 提取细胞核基因名称
    cell_ids = expression["cell_ids"]
    feature_ids = expression["feature_ids"]

    # 创建基础数据，后续添加AnnData对象
    dataset = wrap_data(
        cell_ids=cell_ids,
        feature_ids = feature_ids,
        id=id,
        cell_info=cell_info,
        feature_info=feature_info,
    )
    logger.debug(f"Dataset created: {dataset}")

    # 创建AnnData对象并添加表达矩阵
    dataset = add_expression(dataset, counts["csc"], expression["csc"])
    logger.debug(f"Andata counts and expression added")

    dataset["pydynwrap:with_expression"] = True

    return dataset

def is_wrapper_with_expression(dataset):
    return dataset.get("pydynwrap:with_expression", False)