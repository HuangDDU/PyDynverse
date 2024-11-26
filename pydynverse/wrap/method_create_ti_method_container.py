import os
import docker
from tqdm import tqdm

from .container_get import _container_get_definition
from .method_process_definition import _method_process_definition
from ..util import read_h5, write_h5
from .._logging import logger


def pull_image_with_progress(image_name, tag="latest", logger_func=None):
    # 拉取 Docker 镜像并使用 tqdm 实时显示进度条。
    if logger_func is None:
        # 没有日志函数, 就使用常规print打印
        logger_func = print
    client = docker.from_env()
    try:
        logger_func(f"Try to pull image {image_name}:{tag}...\n")
        api_client = docker.APIClient(
            base_url="unix://var/run/docker.sock")  # 使用 APIClient 获取流式输出
        pull_logs = api_client.pull(
            repository=image_name, tag=tag, stream=True, decode=True)  # 拉取镜像
        progress_bars = {}  # 初始化 tqdm 进度条, 存储每个 layer 的进度条
        for log in pull_logs:
            # 拉取日志为 JSON 格式，解析后展示
            if "status" in log:
                status = log["status"]
                layer_id = log.get("id", None)  # 获取当前的 layer id
                progress_detail = log.get("progressDetail", {})
                current = progress_detail.get("current", 0)  # 当前已完成的字节数
                total = progress_detail.get("total", 0)  # 总字节数
                # 如果有 layer_id 和 total，更新进度条
                if layer_id and total:
                    if layer_id not in progress_bars:
                        # 创建新的进度条
                        progress_bars[layer_id] = tqdm(
                            total=total,
                            desc=f"Layer {layer_id[:12]}",
                            unit="B",
                            unit_scale=True,
                            unit_divisor=1024
                        )
                    progress_bars[layer_id].n = current
                    progress_bars[layer_id].refresh()
                # 如果没有进度信息，打印状态
                elif layer_id:
                    logger_func(f"{status} {layer_id}".strip())
                else:
                    logger_func(f"{status}".strip())
        # 关闭所有进度条
        for bar in progress_bars.values():
            bar.close()
        logger_func(f"Pull {image_name}:{tag} finish")
    except docker.errors.APIError as e:
        logger_func(f"Pull image failed: {e}")
    except Exception as e:
        logger_func(f"Other Error: {e}")
    finally:
        client.close()


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
        # client.images.pull(container_id)
        image_name, tag = container_id.split(":")
        pull_image_with_progress(image_name, tag=tag, logger_func=logger.debug)
        img = client.images.get(container_id)
        logger.debug(f"Docker image({container_id}) loaded")

    # 此时只是有了镜像，容器还没启动

    definition = _container_get_definition(container_id)
    definition["run"] = {"backend": "container", "container_id": container_id}
    definition = _method_process_definition(
        definition=definition, return_function=return_function)
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
    task = inputs  # 表达矩阵，行名列名一起传递
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
    
    log_list = [log.decode("utf-8").strip() for log in container.logs(stream=True)]
    # for log in container.logs(stream=True):
    #     logger.debug(log.decode("utf-8").strip())  # 解码日志并实时打印
    container.wait()  # 当代入口程序执行完成后在执行后续内容
    container.stop()
    container.remove()

    if preproc_meta["verbose"]:
        pass
    
    log = "\n".join(log_list)
    output_h5_filename = f"{tmp_wd}/output.h5"
    if not os.path.exists(output_h5_filename):
        #  没有生成h5文件则docker运行失败, 返回错误信息
        logger.error("Docker Error!!!")
        logger.error(log)
    else:
        logger.debug("Docker Finish")
        logger.debug(log)
        dynverse_docker_output = read_h5(f"{tmp_wd}/output.h5") # 读取输出
        return dynverse_docker_output


def _method_execution_postproc_container(preproc_meta):
    # 容器执行后处理
    pass
