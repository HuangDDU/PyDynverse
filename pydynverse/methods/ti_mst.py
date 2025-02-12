import os

from .method_choose_backend import method_choose_backend
from .function import ti_mst_function


def ti_mst(
    dimred="pca",
    ndim=2
):
    return method_choose_backend(
        package_repository=None,
        package_name=None,
        function_name=None,
        package_version=None,
        definition_filename=f"{os.path.dirname(os.path.abspath(__file__))}/definition/ti_mst_definition.yml",
        run_fun=ti_mst_function,
        container_id="dynverse/ti_mst:v0.9.9.01",
    )(
        dimred=dimred,
        ndim=ndim
    )
