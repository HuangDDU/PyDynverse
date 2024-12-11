from .method_process_definition import _method_load_definition, _method_process_definition


def create_ti_method_py(
    definition,
    run_fun,
    package_required=None,
    package_loaded=None,
    remotes_package=None,
    return_function=True
):
    # 从python脚本创建轨迹推断方法
    definition = _method_load_definition(definition)
    definition["run"] = {
        "backend": "function",
        "run_fun": run_fun,
        "package_required": package_required,
        "package_loaded": package_loaded,
        "remotes_package": remotes_package,
    }

    definition = _method_process_definition(definition=definition, return_function=return_function)

    return definition


def _method_execution_preproc_function(method):
    # 后续函数执行也不使用这些预处理的结果，直接跳过即可
    pass


def _method_execution_execute_function(
        method,
        inputs,
        priors,
        parameters,
        verbose,
        seed,
        preproc_meta
):
    inputs.update({
        "priors": priors,
        "parameters": parameters,
        "verbose": verbose,
        "seed": seed,
    })  # 构造参数
    trajectory = method["run"]["run_fun"](**inputs)  # 执行
    return trajectory


def _method_execution_postproc_function():
    pass
