import pandas as pd

from .wrap_add_trajectory import add_trajectory


def add_cluster_graph(dataset, milestone_network, grouping=None, explicit_splits=False):

    if grouping is None:
        grouping = dataset.get("grouping", None)
    if grouping is None:
        raise "grouping is need for add_cluster_graph"
    if explicit_splits:
        # TODO: 不连通图
        pass

    mn_ft = milestone_network[["from", "to"]]
    both_direction = pd.concat([
        mn_ft.assign(label=mn_ft["from"], percentage=0),
        mn_ft.assign(label=mn_ft["to"], percentage=1)
    ])

    progressions = pd.DataFrame({"cell_id": dataset["cell_ids"], "label": grouping})\
        .merge(both_direction, on="label")\
        .groupby("cell_id")\
        .apply(lambda x: x.sort_values("percentage", ascending=False).iloc[0])\
        .reset_index(drop=True)\
        .drop("label", axis=1)

    trajectory = add_trajectory(
        dataset=dataset,
        milestone_network=milestone_network,
        divergence_regions=None,
        progressions=progressions,
    )
    return trajectory
