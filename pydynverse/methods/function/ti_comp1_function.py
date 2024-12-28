import anndata as ad
import scanpy as sc

from ..._logging import logger
from ...wrap import wrap_data, add_linear_trajectory


def ti_comp1_function(expression, priors, parameters, seed, verbose, cell_ids=None, feature_ids=None, **kwargs):
    # NOTE: 这里还是按照Dynverse的ti_paga的参数，其实这里调用直接传入整个dataset都行
    logger.debug("ti_paga_function executing")
    logger.debug(f"priors: {priors}")
    logger.debug(f"parameters: {parameters}")
    logger.debug(f"seed: {seed}")

    # 1. 数据构造
    adata = ad.AnnData(X=expression)

    # 2. 执行PCA
    sc.pp.pca(adata, n_comps=parameters["ndim"])

    # 3. 结果封装保存
    pseudotime = adata.obsm["X_pca"][:, parameters["component"]-1]
    dataset = wrap_data(cell_ids=cell_ids)
    dataset = add_linear_trajectory(dataset, pseudotime=pseudotime,)
    
    return dataset
