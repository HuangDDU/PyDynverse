from scipy.sparse import csc_matrix
import pandas as pd
from ..util import random_time_string


def wrap_data(
    cell_ids: list,
    feature_ids: list = None,
    id: str = None,
    cell_info: pd.DataFrame = None,
    feature_info: pd.DataFrame = None,
) -> dict:
    if id is None:
        # 生成时间相关的随机字符串作为ID标识数据
        id = random_time_string("data_wrapper")

    assert len(list(set(cell_ids))) == len(
        cell_ids), "cell_ids is not unique"  # 判断cell_ids列表中元素唯一

    if cell_info is None:
        cell_info = pd.DataFrame(index=cell_ids)

    if not feature_ids is None:
        assert len(list(set(feature_ids))) == len(
            feature_ids), "feature_ids is not unique"  # 判断feature_ids列表中元素唯一

        if feature_info is None:
            feature_info = pd.DataFrame(index=feature_ids)

    # AnnData到添加表达矩阵的时候再创建
    # # 创建AnnDatata对象, 方便后续绘图或调用其他方法
    # adata = ad.AnnData(csc_matrix((len(cell_ids), len(feature_ids))))
    # # 添加obs或var的DataFrame
    # if cell_info:
    #     adata.obs = cell_info
    # else:
    #     adata.obs.index = cell_ids
    # if feature_info:
    #     adata.var = feature_info
    # else:
    #     adata.var.index = feature_ids
    # # 标识数据ID
    # adata.uns["id"] = id

    dataset = {
        "id": id,
        "cell_ids": cell_ids,
        "cell_info": cell_info,
        "feature_ids": feature_ids,
        "feature_info": feature_info,
        "pydynwrap:data_wrapper": True
    }

    return dataset


def is_data_wrapper(dataset) -> bool:
    return dataset.get("pydynwrap:data_wrapper", False)