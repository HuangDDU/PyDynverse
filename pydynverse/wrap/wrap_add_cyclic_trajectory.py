import pandas as pd
from .wrap_add_trajectory import add_trajectory
from .wrap_add_pseudotime import process_pseudotime


def add_cyclic_trajectory(
    dataset: dict,
    pseudotime: list,
    directed: bool = False,
    do_scale_minmax: bool = True,
    **kwargs
) -> dict:

    pseudotime = process_pseudotime(dataset, pseudotime)

    if do_scale_minmax:
        pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())
    else:
        assert (pseudotime >= 0).all() and (pseudotime <= 1).all()

    # 构建里程碑网络：A->B, B->C, C->A
    milestone_ids = ["A", "B", "C"]
    milestone_network = pd.DataFrame({
        "from": milestone_ids,
        "to": milestone_ids[1:] + [milestone_ids[0]],
        "length": 1,
        "directed": directed,
        "edge_id": range(len(milestone_ids))
    })

    # 构建细胞过程，3个分段内
    progressions = pd.DataFrame({
        "cell_id": dataset["cell_ids"],
        "time": [3*i for i in pseudotime],
    })
    progressions["edge_id"] = progressions["time"].apply(lambda x: 0 if x <= 1 else 1 if x <= 2 else 2).astype("int")
    progressions = pd.merge(progressions, milestone_network[["from", "to", "edge_id"]], on="edge_id")
    progressions["percentage"] = progressions["time"] - progressions["edge_id"]
    progressions = progressions[["cell_id", "from", "to", "percentage"]].reset_index(drop=True)

    milestone_network = milestone_network[["from", "to", "length", "directed"]]

    trajectory = add_trajectory(
        dataset=dataset,
        milestone_ids=milestone_ids,
        milestone_network=milestone_network,
        divergence_regions=None,
        progressions=progressions,
    )

    return trajectory
