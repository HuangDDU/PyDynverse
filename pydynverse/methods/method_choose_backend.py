from .._settings import settings
from .._logging import logger

from ..wrap import create_ti_method_py, create_ti_method_container
from .method_install_github_tagged_version import install_github_tagged_version
from .method_install_pipy_tagged_version import install_pipy_tagged_version


def method_choose_backend(
        # 这里指定脚本运行的包位置，这其实也是Dynverse自己封装的脚本, 暂时不用这些参数
        package_repository=None,
        package_name=None,
        function_name=None,
        package_version=None,
        # 这里在子包中直接给出
        definition_filename=None,
        run_fun=None,
        # docker镜像位置
        container_id=None,
        backend=None,
):
    # 候选使用wrapper脚本或者Docker容器
    correct_backends = ["function", "container"]

    # 默认从设置中读取
    if backend is None:
        backend = settings["backend"]
    if backend is None:
        # 按照选择后端
        if (not function_name is None) and (container_id is None):
            # 提供了执行函数, 而没有提供docker镜像名
            backend = "function"
        elif (run_fun is None) and (not container_id is None):
            # 没有提供了执行函数, 而提供docker镜像名
            backend = "container"
        elif backend is None:
            # 都没有提供，或者都提供了， 手动选择
            answer = input("""
                    You can run this method as an Python function (1, default) or as a container (2)
                    Which do you want to use?
                    1: Python function [default]
                    2: Container
            """)
            if answer == "2":
                backend = "container"
            else:
                backend = "python_function"
        settings["backend"] = backend
    logger.info(f"backend: {backend}")
    # 执行
    if backend == "python_function":
        # TODO: 下载需要的包, 检查运行环境（暂时不用，直接配置好了环境）
        # install_github_tagged_version() # 从GitHub下载特定版本的包及环境，这里
        # install_pipy_tagged_version(package_name, package_version)  # 从pipy下载特定版本的包
        return create_ti_method_py(definition=definition_filename, run_fun=run_fun)
    elif backend == "container":
        return create_ti_method_container(container_id)  # 创建docker环境
