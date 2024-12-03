import pandas as pd
def convert_milestone_percentages_to_progressions(
    cell_ids,
    milestone_ids,
    milestone_network,
    milestone_percentages
):
    # 第一次连接，直接参考to键相连, 有许多错误cell_id-from对应关系
    df1 = pd.merge(milestone_network, milestone_percentages, left_on="to", right_on="milestone_id")
    # 第二次连接，只保留有效的cell_id-from对应关系
    df2 = pd.merge(df1, milestone_percentages[["cell_id", "milestone_id"]], left_on=["from", "cell_id"], right_on=["milestone_id", "cell_id"])
    progr_part1 = df2[["cell_id", "from", "to", "percentage"]]

    # TODO: progressions第二部分暂时不理解
    progr_part2 = None

    progressions = pd.concat([progr_part1, progr_part2], ignore_index=True)

    return progressions
    
    
def convert_milestone_percentages_to_progressions(
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

