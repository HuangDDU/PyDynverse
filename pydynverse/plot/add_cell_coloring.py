import matplotlib.colors as mcolors
import pandas as pd

from .dummy_proofing import check_groups
from .add_milestone_coloring import add_milestone_coloring
from .mix_colors import mix_colors
from .._logging import logger


def add_cell_coloring(
        cell_positions,
        trajectory,
        color_cells="auto",
        grouping=None,
        groups=None,
        feature_oi=None,
        expression_source="expression",
        pseudotime=None,
        color_milestones=None,
        milestones=None,
        milestone_percentages=None):
    # 选择细胞颜色区分方式
    if color_cells == "auto":
        if grouping:
            logger.info("Coloring by grouping")
            color_cells = "grouping"
        elif feature_oi:
            logger.info("Coloring by expression")
            color_cells = "feature"
        elif (milestones is not None) or (milestone_percentages is not None):
            logger.info("Coloring by milestone")
            color_cells = "milestone"
        elif pseudotime:
            color_cells = "pseudotime"
        else:
            color_cells = "grey"

    result = {}
    if color_cells == "grouping":
        groups = check_groups(grouping, groups)
        color_dict = dict(zip(groups["group_id"], groups["color"]))
        color_scale = [color_dict[i] for i in grouping]
        fill_scale = None
        result["color_dict"] = color_dict
    elif color_cells == "feature":
        # 基因表达值染色
        pass
    elif color_cells == "milestone":
        # 按照grouping的逻辑
        milestone_percentages = trajectory["milestone_percentages"]
        milestone_id_list = trajectory["milestone_ids"]
        if milestones is None:
            milestones = pd.DataFrame([milestone_id_list], columns=["milestone_id"], index=milestone_id_list)
        if "color" not in milestones.columns:
            # 函数传参传入了没有里程碑颜色
            milestones = add_milestone_coloring(milestone_id_list, color_milestones=color_milestones)  # 添加里程碑颜色
        milestone_colors = milestones["color"].apply((lambda x: pd.Series(mcolors.to_rgb(x))))  # 16进制转化为rgb，方便颜色混合
        cell_colors = milestone_percentages.groupby("cell_id").apply(lambda x: mix_colors(x["milestone_id"], x["percentage"], milestone_colors))
        color_scale = cell_colors[cell_positions.index]
        fill_scale = None
    else:
        pass
        # TODO: 其他颜色方式

    result.update({
        "cell_positions": cell_positions,
        "color_scale": color_scale,
        "fill_scale": fill_scale,
        "color_cells": color_cells
    })
    return result
