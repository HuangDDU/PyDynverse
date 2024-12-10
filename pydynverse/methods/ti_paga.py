import os

from .method_choose_backend import method_choose_backend
from .function.ti_paga_function import ti_paga_function


def ti_paga(
    filter_features=True,
    n_neighbors=15,
    n_comps=50,
    n_dcs=15,
    resolution=1,
    embedding_type="fa",
    connectivity_cutoff=0.05
):
    # 选择后端执行, 备选的Docker容器ID, 传参执行
    return method_choose_backend(
        # 这里指定脚本运行的包位置，这其实也是Dynverse自己封装的脚本, 暂时不用这些参数
        package_repository=None,
        package_name=None,
        function_name=None,
        package_version=None,
        # 这里在子包中直接给出
        definition_filename=f"{os.path.dirname(os.path.abspath(__file__))}/definition/ti_paga_definition.yml",
        run_fun=ti_paga_function,
        # docker镜像位置
        container_id="dynverse/ti_paga:v0.9.9.05",
    )(
        filter_features=filter_features,
        n_neighbors=n_neighbors,
        n_comps=n_comps,
        n_dcs=n_dcs,
        resolution=resolution,
        embedding_type=embedding_type,
        connectivity_cutoff=connectivity_cutoff
    )
