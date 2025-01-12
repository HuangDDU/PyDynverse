def group_onto_trajectory_edges(trajectory):
    # lambda y: f"{y['from']}->{y['to']}", axis=1
    def get_trajectory_edges(x):
        x = x.loc[x["percentage"].idxmax()]
        return f"{x['from']}->{x['to']}"
    grouping = trajectory["progressions"].groupby("cell_id").apply(get_trajectory_edges)
    grouping = grouping.loc[trajectory["cell_ids"]].to_list()
    return grouping


def group_onto_nearest_milestones(trajectory):
    def get_nearest_milestone(x):
        return x.loc[x["percentage"].idxmax(), "milestone_id"]
    grouping = trajectory["milestone_percentages"].groupby("cell_id").apply(get_nearest_milestone)
    grouping = grouping.loc[trajectory["cell_ids"]].to_list()
    return grouping
