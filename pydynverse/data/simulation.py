import anndata as ad
import rpy2.robjects as ro
from ..util import rpy2_read  # 导入对应的装饰器， 实现数据自动转换


def load_simulation_data(data_filename="synthetic/dyntoy/bifurcating_1.rds", data_dir="/usr/share/CellFateExplorer/dynbenchmark/data/"):
    #  读取模拟数据，默认读取二分支结构
    r_script = f"""
    dataset <- readRDS("{data_dir}/{data_filename}")
    dataset
    """
    # 这里数据读取有点错乱
    dataset = ro.r(r_script)

    # 添加额外的adata数据
    # adata = ad.AnnData(dataset["counts"])
    # adata.layers["counts"],  adata.layers["expression"] = dataset["counts"], dataset["expression"]
    # 优先添加 expression
    layers = {}
    if "expression" in dataset:
        X = dataset["expression"]
        layers["expression"] = dataset["expression"]
    if "count" in dataset:
        X = dataset["count"]
        layers["count"] = dataset["count"]
    adata = ad.AnnData(X)
    adata.layers = layers
    adata.obs.index = dataset["cell_ids"]

    # 额外操作
    feature_info = dataset.get("feature_info")
    if feature_info is not None:
        if "feature_id" in feature_info.columns:
            dataset["feature_ids"] = dataset["feature_info"]["feature_id"].tolist() # synthetic/dyntoy/bifurcating_1.rds
        if "f" in feature_info.columns:
            dataset["feature_ids"] = dataset["feature_info"]["f"].tolist() # fibroblast-reprogramming_treutlein.rds
        adata.var.index = dataset["feature_ids"]
    dataset["adata"] = adata
    return dataset
