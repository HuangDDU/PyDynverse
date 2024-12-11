import pandas as pd


def convert_progressions_to_milestone_percentages(
    cell_ids,
    milestone_ids,
    milestone_network,
    progressions
):
    # 自环的milstone
    selfs = progressions[progressions["from"] == progressions["to"]]
    selfs = selfs[["cell_id", "from"]].copy().rename(columns={"from": "milestone_id"})
    selfs["percentage"] = 1

    # 非自环的milstone, 后续分别对from和to
    progressions = progressions[~(progressions["from"] == progressions["to"])]

    # 属于from milestone的概率
    froms = progressions[["cell_id", "from", "percentage"]].copy().rename(columns={"from": "milestone_id"})
    froms["percentage"] = 1 - froms["percentage"]

    # 属于to milestone的概率,直接取而对应的列即可
    tos = progressions[["cell_id", "to", "percentage"]].copy().rename(columns={"to": "milestone_id"})

    milestone_percentages = pd.concat([selfs, froms, tos]).reset_index(drop=True)

    return milestone_percentages
