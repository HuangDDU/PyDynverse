import inspect


def inherit_default_params(source_func, param_list=None):
    def decorator(target_func):
        source_params = inspect.signature(source_func).parameters  # 获取源函数的参数信息
        target_params = inspect.signature(
            target_func).parameters  # 获取目标函数的参数信息

        filtered_target_params = {}  # 待过滤的目标函数参数
        if param_list is None:
            # 更新目标函数中的全部参数
            filtered_target_params = target_params
        else:
            # 更新目标函数中的部分参数
            filtered_target_params = {}
            for k, v in target_params.items():
                if k in param_list:
                    filtered_target_params[k] = v

        # 创建一个新的默认参数列表
        new_defaults = []
        for param_name, param in filtered_target_params.items():
            if param_name in source_params:
                # 如果源函数有该参数，并且有默认值，使用源函数的默认值
                if source_params[param_name].default is not inspect.Parameter.empty:
                    new_defaults.append(source_params[param_name].default)
                else:
                    new_defaults.append(
                        param.default if param.default is not inspect.Parameter.empty else None)
            else:
                # 如果源函数没有该参数，使用目标函数的默认值
                new_defaults.append(
                    param.default if param.default is not inspect.Parameter.empty else None)

        # 更新目标函数的 __defaults__
        target_func.__defaults__ = tuple(new_defaults)

        return target_func
    return decorator


if __name__ == "__main__":
    def func2(a=1, b=2, c=3):
        return a + b + c

    @inherit_default_params(func2)
    def func1(b, a):
        return a + b

    @inherit_default_params(func2, ["a"])
    def func3(b, a):
        return a + b

    # 测试
    print(func1())  # 输出: 1(func2.a) + 2(func2.b) = 3
    print(func1(3))  # 输出:  1(func2.a) + 3  = 4

    # print(func3()) # 报错, 此时只有a有参数
    print(func3(b=10))  # 输出: 1(func2.a) + 10 = 11
