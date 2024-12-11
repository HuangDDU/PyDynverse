import pandas as pd
from .wrap_add_trajectory import add_trajectory
import networkx as nx


def add_branch_trajectory(
    dataset,
    branch_network,
    branches,
    branch_progressions,
    **kwargs,
):
    cell_ids = dataset["cell_ids"]
    branch_ids = branches["branch_id"]

    # TODO: 检查branch ids, branch network network和branches，暂时跳过
    branch_network = check_branch_network(branch_ids, branch_network)
    branch = check_branches(branch_ids, branches)
    branch_progressions = check_branch_progressions(cell_ids, branch_ids, branch_progressions)

    # 基于branch_network和branch_progression构建milestone_network和progressions
    milestone_network = pd.DataFrame({
        "from": map(lambda x: f"{x}_from", branch_ids),
        "to": map(lambda x: f"{x}_to", branch_ids),
        "branch_id": branch_ids
    })  # 为每个branch设置起始\终止milestone, 以_from和_to为后缀结尾
    # 构建各个branch的from和to孤立结点, 以及branch交点处的连接图
    milestone_mapper_network = pd.concat(
        [
            # from孤立节点
            pd.DataFrame({
                "from": map(lambda x: f"{x}_from", branch_ids),
                "to": map(lambda x: f"{x}_from", branch_ids),
            }),
            # branch交点处的连接, branch_network A->B, 说明之milestone A_to=B_from, 在一个连通分量里
            pd.DataFrame({
                "from": map(lambda x: f"{x}_to", branch_network["from"]),
                "to": map(lambda x: f"{x}_from", branch_network["to"]),
            }),
            # to孤立节点
            pd.DataFrame({
                "from": map(lambda x: f"{x}_to", branch_ids),
                "to": map(lambda x: f"{x}_to", branch_ids),
            }),
        ]
    )
    # 节点名称转化为连通分量的序号
    mapper = {}
    graph = nx.from_pandas_edgelist(milestone_mapper_network, source="from", target="to")
    connected_components = nx.connected_components(graph)
    for component_index, component in enumerate(connected_components):
        for node in component:
            mapper[node] = str(component_index+1)  # milestone序号从1开始
    milestone_network["from"] = milestone_network["from"].apply(lambda x: mapper[x])
    milestone_network["to"] = milestone_network["to"].apply(lambda x: mapper[x])
    milestone_network = pd.merge(milestone_network, branches, on="branch_id")

    # TODO: 添加额外的自环, 对于PAGA非常重要
    # 对于milstoneA有名称为0的环,则该环扩展为 A -> A-0a -> A-0b -> A
    new_edge_length = milestone_network["length"].mean()/100  # 额外的边长非常小
    # for branch_id in milestone_network[milestone_network["from"]==milestone_network["to"]]["branch_id"]: # 替换这种复杂的写法
    for branch_id in milestone_network.query("`from` == to")["branch_id"]:
        now_branch_index = milestone_network["branch_id"] == branch_id
        now_branch_index = now_branch_index[now_branch_index].index.tolist()
        milestone_id = milestone_network.loc[now_branch_index, "from"].values[0]
        new_milestone_ids = [f"{milestone_id}-{branch_id}{i}" for i in ['a', 'b']]
        milestone_network.loc[now_branch_index, "to"] = new_milestone_ids[0] # A -> A-0a
        mew_milstone_netowrk = pd.DataFrame({
            "from": new_milestone_ids,
            "to": [new_milestone_ids[1], milestone_id],
            "branch_id": "whatever",
            "length": new_edge_length,
            "directed": False
        }) # A-0a -> A-0b -> A
        milestone_network = pd.concat([milestone_network, mew_milstone_netowrk])

    # 直接连接两个branch_progressions与milestone_network, 保留列即可
    progressions = pd.merge(branch_progressions, milestone_network, on="branch_id")[["cell_id", "from", "to", "percentage"]]

    # milestone_network保留列
    milestone_network = milestone_network[["from", "to", "length", "directed"]]

    # 最后调用统一的add_trajectory函数c
    dataset = add_trajectory(
        dataset=dataset,
        milestone_network=milestone_network,
        progressions=progressions
    )
    return dataset


def check_branch_network(branch_ids, branch_network):
    return branch_network


def check_branches(branch_ids, branches):
    return branches


def check_branch_progressions(cell_ids, branch_ids, branch_progressions):
    return branch_progressions
