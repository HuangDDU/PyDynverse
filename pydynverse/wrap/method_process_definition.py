import pandas as pd
from .._logging import logger


def definition(
    method,
    wrapper,
    manuscript=None,
    container=None,
    package=None,
    parameters=None
):
    # 这个函数名字不是很好,容易和之后哪个函数的变量重名
    d = {
        "method": method,
        "wrapper": wrapper,
        "container": container,
        "package": package,
        "manuscript": manuscript,
        "parameters": parameters,
    }

    # 额外添加, 为了方便后续索引, 把参数列表格式改为字典
    # if type(d["parameters"]) is list:
    #     parameter_list = d["parameters"]
    #     parameter_dict = {}
    #     for parameter in parameter_list:
    #         parameter_dict[parameter["id"]] = parameter
    #     d["parameters"] = parameter_dict
    if type(d["parameters"]) is list:
        # 改成DataFrame格式
        df = pd.DataFrame(d["parameters"])
        df.index = df["id"]
        d["parameters"] = df

    # 必选输入和可选输入
    inputs = d["wrapper"]["input_required"]
    inputs = inputs if type(inputs) == list else [inputs]
    if "input_optional" in d["wrapper"]:
        inputs += d["wrapper"]["input_optional"]
    # 参数名称
    params = list(d["parameters"]["id"])

    # 额外的inputs创建, 包括数据和参数
    input_id_list = inputs + params
    required_list = [i in d["wrapper"]["input_required"]
                     for i in input_id_list]
    type_list = []
    for input_id in input_id_list:
        # 输入元素的类型
        if input_id in ["counts", "expression", "expression_future"]:
            type_list.append("expression")
        elif input_id in params:
            type_list.append("parameter")
        else:
            type_list.append("prior_information")
    inputs_df = pd.DataFrame(
        {"input_id": input_id_list, "required": required_list, "type": type_list})
    d["wrapper"]["inputs"] = inputs_df

    return d


def _method_process_definition(definition, return_function):
    if not return_function:
        return definition
    else:

        defaults = get_default_parameters(definition)  # 获取代码函数中的参数默认值

        def param_overrider_fun(**kwargs):
            # 参数覆盖, 接受代码函数中的参数默认值传入参数并覆盖definition.yml
            new_defaults = kwargs
            param_names = list(definition["parameters"].index)
            for param_name, v in new_defaults.items():
                if param_name in param_names:
                    definition["parameters"].loc[param_name, "default"] = v
                else:
                    # 该参数不definition.yml文件里定义
                    logger.error(f"Unknown parameter: {param_name}")
            return definition

        param_overrider_fun.__kwdefaults__ = defaults

        return param_overrider_fun


def convert_definition(definition_raw):
    return definition(
        method=definition_raw["method"],
        wrapper=definition_raw["wrapper"],
        container=definition_raw["container"],
        package=definition_raw["package"] if "package" in definition_raw else None,
        manuscript=definition_raw["manuscript"] if "manuscript" in definition_raw else None,
        parameters=definition_raw["parameters"],
    )


def get_default_parameters(definition):
    # 获取轨迹推断方法的默认参数
    default_parameters = definition["parameters"]["default"]
    return dict(default_parameters)
