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
        "backend": function
    }

    definition = _method_process_definition(
        definition=definition, return_function=return_function)

    return definition


def _method_execution_preproc_function(method):
    run = method["run"]

    # 后续函数调用的执行也不使用


def _method_execution_execute_function(
        method,
        inputs,
        priors,
        parameters,
        verbose,
        seed,
        preproc_meta
):
    args = []
    # 函数方法调用
    trajectory = method["run"]["run_fun"](args)
    return trajectory


def _method_execution_postproc_function():
    pass
