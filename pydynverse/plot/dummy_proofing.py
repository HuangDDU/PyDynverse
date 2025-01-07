import pandas as pd
from .milestone_palette import milestone_palette


def check_milestones(trajectory, milestones, milestone_percentages=None):
    # 没有调用
    milestone_id_list = trajectory["milestone_ids"]
    milestones = pd.DataFrame(data=[milestone_id_list], columns=milestone_id_list)
    return milestones
    


def check_groups(grouping: list, groups: pd.DataFrame):
    """给每个细胞群添加颜色

    Args:
        grouping (list): n_obs
        groups (pd.Dataframe): 已有的群颜色，列名为cell_id, color

    Returns:
        pd.Dataframe: 群颜色，列名为cell_id, color
    """
    if (groups is None) or ("color" not in groups.columns):
        # 没有颜色则手动添加
        group_id = list(set(grouping))  # unique
        color = milestone_palette("auto", len(group_id))  # 自动选择颜色
        groups = {
            "group_id": group_id,
            "color": color
        }
    return groups

