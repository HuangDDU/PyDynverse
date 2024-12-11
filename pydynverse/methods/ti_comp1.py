from .method_choose_backend import method_choose_backend


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
        container_id="dynverse/ti_comp1:v0.9.9.01",
    )(
        dimred=dimred,
        ndim=ndim,
        component=component
    )
