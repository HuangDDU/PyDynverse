import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scvelo as scv

from ..._logging import logger
from ...wrap import wrap_data, add_branch_trajectory


def ti_paga_function(counts, priors, parameters, seed, verbose, cell_ids=None, feature_ids=None, **kwargs):
    # NOTE: 这里还是按照Dynverse的ti_paga的参数，其实这里调用直接传入整个dataset都行
    logger.debug("ti_paga_function executing")
    logger.debug(f"priors: {priors}")
    logger.debug(f"parameters: {parameters}")
    logger.debug(f"seed: {seed}")

    # 1. 数据构造
    start_id = priors["start_id"]
    connectivity_cutoff = parameters["connectivity_cutoff"]
    cluster_key = "clusters"
    adata = ad.AnnData(X=counts)
    adata.obs[cluster_key] = priors["groups_id"]
    if not cell_ids is None:
        adata.obs.index = cell_ids
    if not feature_ids is None:
        adata.var.index = feature_ids

    # 2. 预处理
    scv.pp.filter_and_normalize(adata)
    sc.pp.neighbors(adata, n_neighbors=3)
    sc.tl.diffmap(adata)

    # 3. 方法调用
    # PAGA调用
    sc.tl.paga(adata, groups=cluster_key)
    # 设置起点执行dpt
    adata.uns["iroot"] = np.where(adata.obs.index == start_id)[0][0]
    sc.tl.dpt(adata, n_dcs=2)

    # 4. 结果提取
    # (1)
    epsilon = 1e-3  # 后续对于非常小的数字缩放值
    cell_ids = adata.obs.index.to_list()
    branch_ids = adata.obs[cluster_key].unique().to_list()
    # (2) branches
    branches = pd.DataFrame({
        "branch_id": branch_ids,
        "directed": True,
    })
    branches["length"] = adata.obs[[cluster_key, "dpt_pseudotime"]].groupby(cluster_key).apply(lambda x: x["dpt_pseudotime"].max() - x["dpt_pseudotime"].min()+epsilon).reset_index()[0]
    # (3) branch_network
    branch_network = pd.DataFrame(
        np.triu(adata.uns["paga"]["connectivities"].todense(), k=0),  # 保留上三角矩阵
        index=adata.obs[cluster_key].cat.categories,
        columns=adata.obs[cluster_key].cat.categories
    ).stack().reset_index()
    branch_network.columns = ["from", "to", "length"]
    branch_network = branch_network[branch_network["length"] >= connectivity_cutoff]  # 设置阈值过滤不显著的边
    average_pseudotime_dict = adata.obs.groupby(cluster_key)["dpt_pseudotime"].mean()
    def modify_milestone_network_direction(x):
        if average_pseudotime_dict[x["from"]] <= average_pseudotime_dict[x["to"]]:
            return x
        else:
            x["from"], x["to"] = x["to"], x["from"]
            return x
    branch_network.apply(modify_milestone_network_direction, axis=1) # 调整边的方向
    # 按照from、to的伪时间顺序排列，方便后续milestone编号
    branch_network["from_pseudotime"] = branch_network["from"].apply(lambda x: average_pseudotime_dict[x])
    branch_network["to_pseudotime"] = branch_network["to"].apply(lambda x: average_pseudotime_dict[x])
    branch_network = branch_network.sort_values(["from_pseudotime", "to_pseudotime"])
    branch_network = branch_network[["from", "to"]].reset_index(drop=True)  # 只保留from, to列
    # (4) branch_progressions
    branch_progressions = pd.DataFrame({
        "cell_id": adata.obs.index,
        "branch_id": adata.obs[cluster_key],
        "percentage": adata.obs["dpt_pseudotime"]
    })
    # branch内部按照伪时间排序
    branch_progressions["percentage"] = branch_progressions.groupby("branch_id")["percentage"].apply(lambda x: (x-x.min())/(x.max() - x.min()+epsilon)).values
    branch_progressions

    # 5. 结果封装保存
    dataset = wrap_data(cell_ids=cell_ids)
    dataset = add_branch_trajectory(dataset,
                                    branch_network=branch_network,
                                    branches=branches,
                                    branch_progressions=branch_progressions)
    # dataset.add_timings(checkpoints)

    return dataset
