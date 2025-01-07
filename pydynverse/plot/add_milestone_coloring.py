import pandas as pd
from .milestone_palette import milestone_palette


def add_milestone_coloring(milestone_id_list, color_milestones=None):
    color = milestone_palette("auto", len(milestone_id_list))  # 自动选择颜色
    milestones = pd.DataFrame({
        "milestone_id": milestone_id_list,
        "color": color
    }, index=milestone_id_list)
    return milestones
