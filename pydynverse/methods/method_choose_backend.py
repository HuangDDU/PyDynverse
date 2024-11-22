from ..wrap import create_ti_method_container


def method_choose_backend(
        package_repository,
        package_name,
        function_name,
        package_version,
        container_id,
        backend="container",  # 暂时用docker获得结果
):
    # 候选使用wrapper脚本或者Docker容器
    correct_backends = ["wrapper", "container"]

    # TODO: 手动选择

    # 执行
    if backend == "wrapper":
        pass
        # 下载需要的包
    elif backend == "container":
        return create_ti_method_container(container_id)  # 创建docker环境
