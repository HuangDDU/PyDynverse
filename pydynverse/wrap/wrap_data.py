from scipy.sparse import csc_matrix
import anndata as ad
from ..util import random_time_string


def wrap_data(cell_ids, feature_ids, id=None, cell_info=None, feature_info=None):
    if not id:
        # 生成时间相关的随机字符串作为ID标识数据
        id = random_time_string("data_wrapper")

    # 创建AnnDatata对象
    adata = ad.AnnData(csc_matrix((len(cell_ids), len(feature_ids))))
    # 添加obs或var的DataFrame
    if cell_info:
        adata.obs = cell_info
    else:
        adata.obs.index = cell_ids
    if feature_info:
        adata.var = feature_info
    else:
        adata.var.index = feature_ids
    # 标识数据ID
    adata.uns["id"] = id

    return adata
