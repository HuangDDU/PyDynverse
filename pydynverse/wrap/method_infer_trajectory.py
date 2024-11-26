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
    # 函数执行获得输出, 输出trajectory和summary
    output = _method_execute(
        dataset,
        method,
        parameters,
        give_priors,
        seed,
        verbose,
        return_verbose,
        debug
    )
    design = {
        "model": output["trajectory"],
        "summary": output["summary"],
    }
    return design


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
    # TODO: 除了dataset和method外大量的其他参数，保持与infer_trajectories函数一致
    design = infer_trajectories(
        dataset,
        method,
        parameters=parameters,
        give_priors=give_priors,
        seed=seed,
        verbose=verbose,
        return_verbose=return_verbose,
        debug=debug
    )
    if debug:
        pass
    elif not "model" in design:
        pass
    else:
        return design["model"]
