import inspect

from .wrap_add_expression import get_expression, is_wrapper_with_expression
from .wrap_add_trajectory import is_wrapper_with_trajectory


def add_prior_information(
    dataset,
    start_id=None,
    end_id=None,
    groups_id=None,
    groups_network=None,
    features_id=None,
    groups_n=None,
    start_n=None,
    end_n=None,
    leaves_n=None,
    timecourse_continuous=None,
    timecourse_discrete=None,
    dimred=None,
    verbose=True,
):
    # 提取除了第一个参数dataset和最后一个参数verbose外的所有参数, 比每个手敲优雅, 还容易扩展
    frame = inspect.currentframe()
    var_name_list, _, _, var_dict = inspect.getargvalues(frame)

    prior_information = {}
    for i in range(1, len(var_name_list)-1):
        var_name = var_name_list[i] # 变量名称
        var_value =  var_dict[var_name] # 变量值
        if not var_value is None:
            prior_information[var_name] = var_value

    if is_wrapper_with_trajectory(dataset) and is_wrapper_with_expression(dataset):
        calculated_prior_information = generate_prior_information(
            cell_ids=prior_information["cell_ids"],
            milestone_ids=prior_information["milestone_ids"],
            milestone_network=prior_information["milestone_network"],
            milestone_percentages=prior_information["milestone_percentages"],
            progressions=prior_information["progressions"],
            divergence_regions=prior_information["divergence_regions"],
            expression=get_expression(dataset),
            feature_info=prior_information["feature_info"],
        )
        prior_information.update(
            calculated_prior_information)  # 计算得到的额外的先验知识更新

    dataset["prior_information"] = prior_information

    return dataset


def generate_prior_information(
    cell_ids,
    milestone_ids,
    milestone_network,
    milestone_percentages,
    progressions,
    divergence_regions,
    expression,
    feature_info=None,
    cell_info=None,
    marker_fdr=0.005,
    given=None,
    verbose=False
):
    return {}
