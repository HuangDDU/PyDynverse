import docker

from .container_get import _container_get_definition
from .method_process_definition import _method_process_definition
from ..util import read_h5, write_h5
from .._logging import logger


def create_ti_method_container(
    container_id,
    # pull_if_needed=True, # 基本不用，没有镜像肯定就得重新拉取
    return_function=True
):
    # 检查镜像

    # 暂时只支持docker容器，不支持singularity容器, 跳过config的选择
    # docker容器版本查看
    client = docker.from_env()

    # 检查镜像是否存在
    try:
        img = client.images.get(container_id)
        logger.debug(f"Docker image({container_id}) loaded")
    except Exception as e:
        # 镜像不存在，需要使用代码下载
        logger.debug(e)
        logger.debug(f"Docker image({container_id}) was not found")
        logger.debug(f"Try to pull docker image")
        client.images.pull(container_id)
        img = client.images.get(container_id)
        logger.debug(f"Docker image({container_id}) loaded")

    # 此时只是有了镜像，容器还没启动

    definition = _container_get_definition(container_id)
    definition["run"] = {"backend": "container", "container_id": container_id}
    definition = _method_process_definition(definition = definition, return_function = return_function)
    return definition


def _method_execution_preproc_container(
        method,
        inputs,
        tmp_wd,
        priors,
        parameters,
        verbose,
        seed,
        debug
):
    # 容器执行前处理， 需要的数据h5数据
    # 构造R对象的
    task = inputs # 表达矩阵，行名列名一起传递
    # task["cell_ids"] = inputs["cell_ids"]  # 表达矩阵的行名
    # task["feature_ids"] = inputs["feature_ids"]  # 表达矩阵的列明
    task["priors"] = priors
    task["parameters"] = parameters
    task["verbose"] = verbose
    task["seed"] = seed

    write_h5(task, f"{tmp_wd}/input.h5")

    path = {}
    path["prior_names"] = "1"
    path["debug"] = debug
    path["verbose"] = verbose

    return path


def _method_execution_execute_container(method, preproc_meta, tmp_wd):
    # 容器执行

    # 构建执行的参数
    args = ["--dataset", "/ti/input.h5", "--output", "/ti/output.h5"]
    # 是否开启调试模式
    if preproc_meta["debug"]:
        args.append("--debug")
    # 是否使用先验知识, 这里暂时不启用，否则容器无法读取
    # if preproc_meta["prior_names"]:
    #     args+=["--use_priors", "all"]
    logger.debug(f"docker args: {args}")

    # 启动容器
    client = docker.from_env()
    container = client.containers.run(
        image=method["run"]["container_id"],
        command=args,  # 入口点程序参数
        volumes=[f"{tmp_wd}:/ti"],
        working_dir="/ti/workspace",
        detach=True,
    )
    logger.debug(container.logs())
    container.wait()  # 当代入口程序执行完成后在执行后续内容
    container.stop()
    container.remove()
    logger.debug("Docker Finish")

    if preproc_meta["verbose"]:
        pass

    # 读取输出
    dynverse_docker_output = read_h5(f"{tmp_wd}/output.h5")

    return dynverse_docker_output


def _method_execution_postproc_container(preproc_meta):
    # 容器执行后处理
    pass
