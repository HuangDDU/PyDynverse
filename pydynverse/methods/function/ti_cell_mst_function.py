import numpy as np
import pandas as pd
import networkx as nx
import anndata as ad
import scanpy as sc

from ..._logging import logger
from ...wrap import wrap_data, add_cell_graph


def ti_cell_mst_function(counts, priors, parameters, seed, verbose, cell_ids=None, feature_ids=None, **kwargs):
    # NOTE: 这里还是按照Dynverse的ti_paga的参数，其实这里调用直接传入整个dataset都行
    logger.debug("ti_paga_function executing")
    logger.debug(f"priors: {priors}")
    logger.debug(f"parameters: {parameters}")
    logger.debug(f"seed: {seed}")

    # 1. 数据构造
    adata = ad.AnnData(X=counts)
    adata.obs.reset_index(drop=True, inplace=True)

    # 2. 执行PCA
    sc.pp.pca(adata, n_comps=parameters["ndim"])

    # 3. 直接构建细胞间的最小生成树
    sc.pp.neighbors(adata)
    G = nx.from_scipy_sparse_array(adata.obsp["distances"])  # 从稀疏矩阵构造图
    cell_mst = nx.minimum_spanning_tree(G, weight="weight")

    # 4. 提取并封装结果
    cell_graph = nx.to_pandas_edgelist(cell_mst, source="from", target="to").rename(columns={"weight": "length"})
    cell_graph["from"] = cell_graph["from"].apply(lambda x: cell_ids[x])
    cell_graph["to"] = cell_graph["to"].apply(lambda x: cell_ids[x])
    to_keep = pd.Series(data=True, index=cell_ids)

    dataset = wrap_data(cell_ids=cell_ids)

    dataset = add_cell_graph(
        dataset=dataset,
        cell_graph=cell_graph,
        to_keep=to_keep
    )

    dataset["adata"] = adata

    return dataset
