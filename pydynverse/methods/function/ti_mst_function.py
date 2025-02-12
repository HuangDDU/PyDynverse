import numpy as np
import pandas as pd
import networkx as nx
import anndata as ad
import scanpy as sc

from sklearn.metrics.pairwise import pairwise_distances

from ..._logging import logger
from ...wrap import wrap_data, add_dimred_projection


def ti_mst_function(expression, priors, parameters, seed, verbose, cell_ids=None, feature_ids=None, **kwargs):
    # NOTE: 这里还是按照Dynverse的ti_paga的参数，其实这里调用直接传入整个dataset都行
    logger.debug("ti_paga_function executing")
    logger.debug(f"priors: {priors}")
    logger.debug(f"parameters: {parameters}")
    logger.debug(f"seed: {seed}")

    # 1. 数据构造
    adata = ad.AnnData(X=expression)
    adata.obs.reset_index(drop=True, inplace=True)

    # 2. 执行PCA
    sc.pp.pca(adata, n_comps=parameters["ndim"])
    X_emb = adata.obsm["X_pca"]
    # n_comps = parameters.get("n_comps", 10)
    # n_gene = expression.shape[1]
    # if n_gene < n_comps:
    #     # 如果基因数小于n_comps，不降维
    #     n_comps = n_gene
    # else:
    #     sc.pp.pca(adata, n_comps=parameters["ndim"])

    # 3. 聚类细胞，中心点作为里程碑
    # （1）为了方便，这里直接调用scanpy的聚类方法leiden
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)
    # （2）计算聚类中心的低维坐标
    centers = np.array(list(adata.obs.groupby("leiden").apply(lambda x: X_emb[list(x.index)].mean(axis=0))))
    milestone_ids = [f"M{i}"for i in range(centers.shape[0])]
    centers = pd.DataFrame(centers, index=milestone_ids)
    # （3）计算聚类中心间的距离
    dis = pd.DataFrame(pairwise_distances(centers, metric="euclidean"), index=milestone_ids, columns=milestone_ids)
    disdf = pd.DataFrame(data=dis.unstack().reset_index().values, columns=["from", "to", "weight"])  # 转化为长数据

    # 4. 里程碑之间计算距离并构建最小生成树作为里程碑网络
    G = nx.from_pandas_edgelist(disdf, source="from", target="to", edge_attr="weight")
    mst = nx.minimum_spanning_tree(G, weight="weight")
    milestone_network = nx.to_pandas_edgelist(mst)
    milestone_network.rename(columns={"source": "from", "target": "to", "weight": "length"}, inplace=True)
    milestone_network["directed"] = False

    # 5. 结果封装保存
    dataset = wrap_data(cell_ids=cell_ids)

    comp_ids = [f"comp_{i+1}"for i in range(centers.shape[1])]
    dimred = pd.DataFrame(X_emb, index=cell_ids, columns=comp_ids)
    centers.columns = comp_ids
    dimred_milestones = centers
    
    dataset = add_dimred_projection(
        dataset=dataset,
        milestone_network=milestone_network,
        dimred=dimred,
        dimred_milestones=dimred_milestones,
    )
    dataset["adata"] = adata

    return dataset
