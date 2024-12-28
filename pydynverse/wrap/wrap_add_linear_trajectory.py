import pandas as pd
from .._logging import logger

from .wrap_add_trajectory import add_trajectory
from .wrap_add_pseudotime import process_pseudotime


def add_linear_trajectory(
    dataset: dict,
    pseudotime: list,
    directed: bool = False,
    do_scale_minmax: bool = True,
    **kwargs
) -> dict:

    pseudotime = process_pseudotime(dataset, pseudotime)

    # 最大最小归一化到[0, 1], 不用调用函数，直接一行写完
    if do_scale_minmax:
        pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())
    else:
        assert (pseudotime >= 0).all() and (pseudotime <= 1).all()

    milestone_ids = ["milestone_begin", "milestone_end"]
    # milestone_network datframe构建，长度简单为1, 一般为有向图
    milestone_network = pd.DataFrame({
        "from": milestone_ids[0],
        "to": milestone_ids[1],
        "length": 1,
        "directed": directed,
    }, index=[0]) # 全部传入标量，需要index=[0]表示只有一行来指定样本数量
    # progressions dataframe构建， percentage即为pseudotime
    progressions = pd.DataFrame({
        "cell_id": dataset["cell_ids"],
        "from": milestone_ids[0],
        "to": milestone_ids[1],
        "percentage": pseudotime,
    })

    # 转化为direct添加轨迹推断结果
    dataset = add_trajectory(
        dataset=dataset,
        milestone_ids=milestone_ids,
        milestone_network=milestone_network,
        divergence_regions=None,
        progressions=progressions,
        pseudotime=pseudotime,
        **kwargs
    )
    return dataset
