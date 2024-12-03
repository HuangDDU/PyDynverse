from .wrap_add_trajectory import add_trajectory

def add_branch_trajectory(
    dataset,
    branch_network,
    branches,
    branch_progressions,
):
    cell_ids = dataset["cell_ids"]
    branch_ids = dataset["branch_ids"]

    # 检查branch ids, branch network network和branches
    branch_network = check_branch_network()
    branch = check_branch()
    branch_progressions = check_branch_progressions()

    # 基于branch_network和branch_progression构建milestone_network和progressions
    milestone_network = branch_network
    progressions = branch_progressions

    # 最后调用统一的add_trajectory函数c
    add_trajectory(
        dataset=dataset,
        milestone_network = milestone_network,
        progressions = progressions
    )


def check_branch_network():
    pass


def check_branch():
    pass


def check_branch_progressions():
    pass
