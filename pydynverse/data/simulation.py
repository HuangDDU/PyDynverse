import anndata as ad
import rpy2.robjects as ro
from ..util import rpy2_read  # 导入对应的装饰器， 实现数据自动转换


def load_simulation_data(data_filename="synthetic/dyntoy/bifurcating_1.rds", data_dir="/home/huang/RCode/scrna_tools/dynbenchmark/data"):
    #  读取模拟数据，默认读取二分支结构
    r_script = f"""
    dataset <- readRDS("{data_dir}/{data_filename}")
    dataset
    """
    # 这里数据读取有点错乱
    dataset = ro.r(r_script)

    # 额外操作
    dataset["feature_ids"] = dataset["feature_info"]["feature_id"].tolist()
    # 添加额外的adata数据
    adata = ad.AnnData(dataset["counts"])
    adata.layers["counts"],  adata.layers["expression"] = dataset["counts"], dataset["expression"]
    adata.obs.index = dataset["cell_ids"]
    adata.var.index = dataset["feature_ids"]
    dataset["adata"] = adata
    return dataset
