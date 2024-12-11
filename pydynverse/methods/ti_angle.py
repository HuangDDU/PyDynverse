from .method_choose_backend import method_choose_backend


def ti_angle(
        dimred="pca"
):
    return method_choose_backend(
        package_repository=None,
        package_name=None,
        function_name=None,
        package_version=None,
        container_id="dynverse/ti_angle:v0.9.9.02",
    )(
        dimred=dimred
    )
