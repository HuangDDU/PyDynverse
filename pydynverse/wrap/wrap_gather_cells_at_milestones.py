from .wrap_add_trajectory import add_trajectory


def gather_cells_at_milestones(trajectory):
    # 把细胞聚合到最近的里程碑上
    milestone_percentages = trajectory["milestone_percentages"]
    
    # 选择细胞所属概率最大的milstone保留且概率为1
    max_percentage_index = milestone_percentages.groupby("cell_id")["percentage"].idxmax() # 把细胞聚合到最近的里程碑上
    milestone_percentages_max = milestone_percentages.loc[max_percentage_index] # 选择这些行
    milestone_percentages_max["percentage"] = 1 # 将这些行的 percentage 列值设置为 1

    trajectory = add_trajectory(
        trajectory,
        milestone_network=trajectory["milestone_network"],
        divergence_regions=trajectory["divergence_regions"],
        milestone_percentages=milestone_percentages_max,
        trajectory_cell_positioning="milestones"
    )
    return trajectory
