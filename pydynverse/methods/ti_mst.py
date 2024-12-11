from .method_choose_backend import method_choose_backend


def ti_mst(
    dimred="pca",
    ndim=2
):
    return method_choose_backend(
        package_repository=None,
        package_name=None,
        function_name=None,
        package_version=None,
        container_id="dynverse/ti_mst:v0.9.9.01",
    )(
        dimred=dimred,
        ndim=ndim
    )
