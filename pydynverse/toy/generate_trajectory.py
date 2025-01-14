import pandas as pd
from ..wrap import wrap_data, add_trajectory


def generate_trajectory(milestone_network: pd.DataFrame, id: str = ""):
    # TODO: 完整实现需要参考dyntoy.generate_trajectory
    if "length" not in milestone_network.columns:
        milestone_network["length"] = 1
    if "directed" not in milestone_network.columns:
        milestone_network["directed"] = True

    # 简单给第一个milestone上添加1个细胞即可
    cell_ids = ["a"]
    dataset = wrap_data(cell_ids, id=id)
    divergence_regions = pd.DataFrame(columns=["divergence_id", "milestone_id", "is_star"])
    milestone_percentages = pd.DataFrame(
        data=[["a", milestone_network.iloc[0, 0], 1.0],],
        columns=["cell_id", "milestone_id", "percentage"])

    trajectory = add_trajectory(
        dataset=dataset,
        milestone_network=milestone_network,
        divergence_regions=divergence_regions,
        milestone_percentages=milestone_percentages
    )

    return trajectory
