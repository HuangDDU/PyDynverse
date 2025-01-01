from .method_execute import _method_execute

from ..util import inherit_default_params


def infer_trajectories(dataset,
                       method,
                       parameters=None,
                       give_priors=None,
                       seed=None,
                       verbose=None,
                       return_verbose=None,
                       debug=None):
    # 这里认为可以同时执行多个dataset和method
    if not isinstance(dataset, list):
        dataset = [dataset]
    if not isinstance(method, list):
        method = [method]

    output = []
    for i in range(len(dataset)):
        for j in range(len(method)):
            tmp_output = _method_execute(
                dataset[i],
                method[j],
                parameters,
                give_priors,
                seed,
                verbose,
                return_verbose,
                debug
            )
            output.append({
                "dataset_ix": i,
                "method_ix": j,
                "dataset_id": dataset[i]["id"],
                "method_id": method[j]["method"]["id"],
                "method_name": method[j]["method"]["name"],
                "model": tmp_output["trajectory"],
                "summary": tmp_output["summary"],
            })

    return output


@inherit_default_params(infer_trajectories)
def infer_trajectory(
    dataset,
    method,
    parameters,
    give_priors,
    seed,
    verbose,
    return_verbose,
    debug

):
    #  除了dataset和method外大量的其他参数，保持与infer_trajectories函数一致
    output = infer_trajectories(
        dataset,
        method,
        parameters=parameters,
        give_priors=give_priors,
        seed=seed,
        verbose=verbose,
        return_verbose=return_verbose,
        debug=debug
    )
    # 这里获得的是单个方法上的单个轨迹推断结果
    if debug:
        pass
    elif not "model" in output[0]:
        pass
    else:
        return output[0]["model"]
