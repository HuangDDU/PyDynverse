from .._logging import logger

from .dummy_proofing import check_groups

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
    if  color_cells == "auto":
        if grouping:
            logger.info("Coloring by grouping")
            color_cells = "grouping"
        elif feature_oi:
            logger.info("Coloring by expression")
            color_cells = "feature"
        elif milestones or milestone_percentages:
            logger.info("Coloring by milestone")
            color_cells <- "milestone"
        elif pseudotime:
            color_cells <- "pseudotime"
        else:
            color_cells <- "grep"

    result = {}
    if color_cells == "grouping":
        groups = check_groups(grouping, groups)
        color_dict = dict(zip(groups["group_id"], groups["color"]))
        color_scale = [color_dict[i] for i in grouping]
        fill_scale = None
        result["color_dict"] = color_dict
    else:
        pass

    result.update({
        "cell_positions": cell_positions,
        "color_scale": color_scale,
        "fill_scale": fill_scale,
        "color_cells": color_cells
    })
    return result