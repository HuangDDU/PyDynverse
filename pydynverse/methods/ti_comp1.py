import os

from .method_choose_backend import method_choose_backend
from .function import ti_comp1_function


def ti_comp1(
    dimred="pca",
    ndim=2,
    component=1
):
    # 选择后端执行, 备选的Docker容器ID, 传参执行
    return method_choose_backend(
        package_repository=None,
        package_name=None,
        function_name=None,
        package_version=None,
        # 这里在子包中直接给出
        definition_filename=f"{os.path.dirname(os.path.abspath(__file__))}/definition/ti_comp1_definition.yml",
        run_fun=ti_comp1_function,
        container_id="dynverse/ti_comp1:v0.9.9.01",
    )(
        dimred=dimred,
        ndim=ndim,
        component=component
    )
