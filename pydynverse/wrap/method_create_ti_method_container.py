import docker

from .._logging import logger

def create_ti_method_container(
        container_id,
        # pull_if_needed=True, # 基本不用，没有镜像肯定就得重新拉取
        return_function=True
    ):
    
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
    return img
