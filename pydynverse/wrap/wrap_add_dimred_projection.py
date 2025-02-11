from ..util import project_to_segments
from .wrap_add_trajectory import add_trajectory
from .wrap_add_dimred import add_dimred


def add_dimred_projection(dataset, milestone_network, dimred, dimred_milestones, grouping=None):
    if grouping is None:
        proj = project_to_segments(
            x=dimred,
            segment_start=dimred_milestones.loc[milestone_network["from"],],
            segment_end=dimred_milestones.loc[milestone_network["to"],],
        )
        progressions = None
    else:
        progressions = None

    dimred_segment_progressions = None
    dimred_segment_points = None

    out = dataset
    # out = add_trajectory(
    #     dataset=dataset,
    #     milestone_network=milestone_network,
    #     divergence_regions=None,
    #     progressions=progressions,
    # )

    # out = add_dimred(
    #     dataset=out,
    #     dimred=dimred,
    #     dimred_milestones=dimred_milestones
    # )

    return out
