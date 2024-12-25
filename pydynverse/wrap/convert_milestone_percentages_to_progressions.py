import pandas as pd


def convert_milestone_percentages_to_progressions(
    cell_ids,
    milestone_ids,
    milestone_network,
    milestone_percentages
):
    # 第一部分，有两个milestone或者更多milestone的细胞
    # part1: for cells that have 2 or more milestones
    # 第一次连接，直接参考to键相连, 有许多错误cell_id-from对应关系
    # first merge based on "to" key result in many invalid cell_id-form relationship
    df1 = pd.merge(milestone_network, milestone_percentages, left_on="to", right_on="milestone_id")
    # 第二次连接，只保留有效的cell_id-from对应关系
    # second merge based on "to" key
    df2 = pd.merge(df1, milestone_percentages[["cell_id", "milestone_id"]], left_on=["from", "cell_id"], right_on=["milestone_id", "cell_id"])
    # TODO: 检测是否可以把这两部合并操作同时做
    # TODO: check if the two step merge can be done simutaneously
    progr_part1 = df2[["cell_id", "from", "to", "percentage"]]

    
    # 第二部分，只有一个milestone的细胞
    # TODO: 暂时只是简单的保留只有一个milestone的cell
    progr_part2 = milestone_percentages.groupby("cell_id").filter(lambda x: len(x) == 1)
    progr_part2["from"] = progr_part2["milestone_id"]
    progr_part2["to"] = progr_part2["milestone_id"]
    progr_part2 = progr_part2[["cell_id", "from", "to", "percentage"]]

    # progressions = pd.concat([progr_part1], ignore_index=True)
    progressions = pd.concat([progr_part1, progr_part2], ignore_index=True).reset_index(drop=True)

    return progressions


def convert_milestone_percentages_to_progressions_me(
    cell_ids,
    milestone_ids,
    milestone_network,
    milestone_percentages
):
    # 确定细胞的边集合（两个元素的普通边，多个元素的延迟承诺区域）
    cell_edge = pd.DataFrame(milestone_percentages.groupby("cell_id")["milestone_id"].agg(set))
    cell_edge.columns = ["edge"]
    cell_edge["cell_id"] = cell_edge.index

    # 确定细胞所属的边（只有两个元素的普通边）
    progressions = []
    for cell_index, cell_row in cell_edge.iterrows():
        for milestone_index, milestone_row in milestone_network.iterrows():
            if {milestone_row["from"], milestone_row["to"]}.issubset(cell_row["edge"]):
                progressions.append([cell_row["cell_id"], milestone_row["from"], milestone_row["to"]])
    progressions = pd.DataFrame(progressions, columns=["cell_id", "from", "to"])

    # 对应链接获得progression
    progressions = pd.merge(progressions, milestone_percentages, left_on=["cell_id", "to"], right_on=["cell_id", "milestone_id"])
    progressions = progressions[["cell_id", "from", "to", "percentage"]]

    return progressions
