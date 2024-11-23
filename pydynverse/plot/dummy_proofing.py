from .milestone_palette import milestone_palette

def check_milestones(trajectory, milestones, milestone_percentages = None):
    pass

def check_groups(grouping, groups):
    group_id = list(set(grouping)) # unique
    color = milestone_palette("auto", len(group_id)) # 自动选择颜色
    groups = {
        "group_id" : group_id,
        "color": color
    }
    return groups