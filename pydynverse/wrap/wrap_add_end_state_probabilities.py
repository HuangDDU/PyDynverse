import pandas as pd
from .wrap_add_linear_trajectory import add_linear_trajectory
from .wrap_add_trajectory import add_trajectory


def add_end_state_probabilities(
        dataset: dict,
        end_state_probabilities: pd.DataFrame,
        pseudotime: list,
        do_scale_minmax: bool = True):

    if do_scale_minmax:
        pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())
    if end_state_probabilities.shape[1] == 1:
        # 只有一个终端状态，就是线性轨迹了
        trajectory = add_linear_trajectory(
            dataset=dataset,
            pseudotime=pseudotime,
            directed=True,
            do_scale_minmax=do_scale_minmax,
        )
    else:
        # 多个终端状态， 构建里程碑网络
        start_milestone_id = "milestone_begin"  # 起始点是一个完全虚拟点
        end_milestone_ids = end_state_probabilities.columns[1:].tolist()  # 终端点从列名中提取, 默认第一列为cell_id
        milestone_ids = [start_milestone_id] + end_milestone_ids

        # 起始点作为中心的星型里程碑网络
        milestone_network = pd.DataFrame({
            "from": start_milestone_id,
            "to": end_milestone_ids,
            "length": 1,
            "directed": True
        })

        # 添加发散区域，由所有里程碑节点共同构成构成
        divergence_regions = pd.DataFrame({
            "milestone_id": milestone_ids,
            "divergence_id": "D",
            "is_start": pd.Series(milestone_ids) == start_milestone_id
        })

        pseudotime = pd.Series(pseudotime, index=end_state_probabilities["cell_id"])
        progressions = end_state_probabilities.melt(id_vars=["cell_id"], var_name="to", value_name="percentage")
        progressions["from"] = start_milestone_id
        progressions["percentage"] = progressions.groupby("cell_id")["percentage"].transform(lambda x: x / x.sum() * pseudotime[x.name])  # 缩放使其之和为1，暂时不理解这个
        progressions = progressions[["cell_id", "from", "to", "percentage"]]

        trajectory = add_trajectory(
            dataset=dataset,
            milestone_network=milestone_network,
            divergence_regions=divergence_regions,
            progressions=progressions
        )

    return trajectory
