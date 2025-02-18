import os

from .method_choose_backend import method_choose_backend
from .function import ti_cell_mst_function


def ti_cell_mst(
    dimred="pca",
    ndim=2
):
    return method_choose_backend(
        package_repository=None,
        package_name=None,
        function_name=None,
        package_version=None,
        definition_filename=f"{os.path.dirname(os.path.abspath(__file__))}/definition/ti_cell_mst_definition.yml",
        run_fun=ti_cell_mst_function,
        container_id="***",
        backend="python_function" # 自己写的baseline，只能用python脚本调用
    )(
        dimred=dimred,
        ndim=ndim
    )
