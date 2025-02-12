import pandas as pd

from ..util import project_to_segments
from .wrap_add_trajectory import add_trajectory
from .wrap_add_dimred import add_dimred


def add_dimred_projection(dataset, milestone_network, dimred, dimred_milestones, grouping=None):

    cell_ids = dataset["cell_ids"]

    if grouping is None:
        # 聚类标签没有给定，直接投影细胞到线段上

        # 计算投影proj字典：
        #   x_proj: 投影点的坐标
        #   distance: 源点到投影点距离
        #   segement: 投影到第几条边上
        #   progression: 投影点在对应边上的比例
        proj = project_to_segments(
            x=dimred,
            segment_start=dimred_milestones.loc[milestone_network["from"],],
            segment_end=dimred_milestones.loc[milestone_network["to"],],
        )

        # 投影的progression保存
        progressions = milestone_network.iloc[proj["segment"]-1][["from", "to"]]
        progressions["cell_id"] = cell_ids
        progressions["percentage"] = proj["progression"]
        progressions = progressions[["cell_id", "from", "to", "percentage"]].reset_index(drop=True)
    else:
        # TODO: 给定了聚类标签，把细胞投影到所属聚类对应的线段上
        progressions = None

    # 降维相关的保存
    def expand_row(row):
        new_df = pd.concat([row]*2, axis=1).T
        new_df["percentage"] = [0, 1]
        new_df["milestone_id"] = [row["from"], row["to"]]
        return new_df

    dimred_segment_progressions = pd.concat(list(milestone_network[["from", "to"]].apply(lambda x: expand_row(x), axis=1))).reset_index(drop=True)
    dimred_segment_points = dimred_milestones.loc[dimred_segment_progressions["milestone_id"]]
    dimred_segment_progressions = dimred_segment_progressions[["from", "to", "percentage"]]

    out = add_trajectory(
        dataset=dataset,
        milestone_network=milestone_network,
        divergence_regions=None,
        progressions=progressions,
    )

    out = add_dimred(
        dataset=out,
        dimred=dimred,
        dimred_milestones=dimred_milestones,
        dimred_segment_points=dimred_segment_points,
        dimred_segment_progressions=dimred_segment_progressions
    )

    if grouping is not None:
        out["grouping"] = grouping

    return out
