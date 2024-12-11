import tempfile
import yaml
import docker

from .method_process_definition import convert_definition

# NOTE: 这里其实放在method_process_definition更好
def _container_get_definition(container_id):
    # 提取docker镜像的definition.yml文件, 包括对于方法的描述/参数等
    with tempfile.TemporaryDirectory() as tmp_wd:
        # 启动容器
        client = docker.from_env()
        container = client.containers.run(
            entrypoint="cp /code/definition.yml /copy_mount/definition.yml",  # 复制挂载文件夹/copy_mount
            image=container_id,
            volumes=[f"{tmp_wd}:/copy_mount"],
            detach=True,
        )
        container.wait()
        container.stop()
        container.remove()
        # 读取yml文件并解析
        with open(f"{tmp_wd}/definition.yml", 'r') as file:
            definition_raw = yaml.safe_load(file)

    definition = convert_definition(definition_raw)

    # # 为了方便后续索引, 把参数列表格式改为字典c
    # parameter_list = definition["parameters"]
    # parameter_dict = {}
    # for parameter in parameter_list:
    #     parameter_dict[parameter["id"]] = parameter
    # definition["parameters"] = parameter_dict
    return definition
