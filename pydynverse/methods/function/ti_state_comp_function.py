import numpy as np
import pandas as pd
import networkx as nx
import anndata as ad
import scanpy as sc

from ..._logging import logger
from ...wrap import wrap_data, add_end_state_probabilities


def ti_state_comp_function(counts, priors, parameters, seed, verbose, cell_ids=None, feature_ids=None, **kwargs):
    # NOTE: 这里还是按照Dynverse的ti_paga的参数，其实这里调用直接传入整个dataset都行
    logger.debug("ti_paga_function executing")
    logger.debug(f"priors: {priors}")
    logger.debug(f"parameters: {parameters}")
    logger.debug(f"seed: {seed}")

    # 1. 数据构造
    adata = ad.AnnData(X=counts)
    adata.obs.reset_index(drop=True, inplace=True)

    # 2. 执行PCA
    ndim = parameters["ndim"]
    sc.pp.pca(adata, n_comps=ndim)

    # 3. 结果封装保存
    pseudotime = adata.obsm["X_pca"][:, parameters["component"]-1]
    comp_column_list = [f"comp_{i}" for i in range(1, ndim+1)]
    end_state_probabilities = pd.DataFrame(
        columns=comp_column_list,
        data=adata.obsm["X_pca"],
        index=cell_ids,
    )
    end_state_probabilities["cell_id"] = cell_ids
    end_state_probabilities = end_state_probabilities[["cell_id"] + comp_column_list]

    dataset = wrap_data(cell_ids=cell_ids)
    dataset = add_end_state_probabilities(
        dataset,
        end_state_probabilities=end_state_probabilities,
        pseudotime=pseudotime,
    )

    return dataset
