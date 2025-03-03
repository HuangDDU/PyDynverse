import os

from .method_choose_backend import method_choose_backend
from .function import ti_angle_function


def ti_angle(
        dimred="pca"
):
    return method_choose_backend(
        package_repository=None,
        package_name=None,
        function_name=None,
        package_version=None,
        definition_filename=f"{os.path.dirname(os.path.abspath(__file__))}/definition/ti_angle_definition.yml",
        run_fun=ti_angle_function,
        container_id="dynverse/ti_angle:v0.9.9.02",
    )(
        dimred=dimred
    )
