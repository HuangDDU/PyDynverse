from ..wrap import create_ti_method_container
from .method_install_github_tagged_version import install_github_tagged_version
from .method_install_pipy_tagged_version import install_pipy_tagged_version


def method_choose_backend(
        package_repository,
        package_name,
        function_name,
        package_version,
        container_id,
        backend="container",  # 暂时用docker获得结果
):
    # 候选使用wrapper脚本或者Docker容器
    correct_backends = ["python_wrapper", "container"]

    
    if (not (function_name is None)) and (container_id is None):
        # 提供了方法名, 而没有提供docker镜像名
        backend = "python_wrapper"
    elif (function_name is None) and (not (container_id is None)):
        # 没有提供了方法名, 而提供docker镜像名
        backend = "container"
    elif backend is None:
        # TODO: 手动选择
        answer = "2"
        if answer == "2":
            backend = "container"
        else:
            backend = "python_wrapper"
    else:
        # 否则，默认选择wrapper脚本执行
        backend = "python_wrapper"

    backend = "container"
    # 执行
    if backend == "python_wrapper":
        # install_github_tagged_version() # 从GitHub下载特定版本的包及环境，这里
        install_pipy_tagged_version(package_name, package_version) # 从pipy下载特定版本的包
        # 下载需要的包, 检查运行环境
    elif backend == "container":
        return create_ti_method_container(container_id)  # 创建docker环境
