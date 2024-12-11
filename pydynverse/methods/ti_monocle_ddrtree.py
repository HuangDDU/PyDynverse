from .method_choose_backend import method_choose_backend


def ti_monocle_ddrtree(
    reduction_method="DDRTree",
    max_components=2,
    norm_method="log",
    auto_param_selection=True,
    filter_features=True,
    filter_features_mean_expression=0.1
):
    # 选择后端执行, 备选的Docker容器ID, 传参执行
    return method_choose_backend(
        package_repository=None,
        package_name=None,
        function_name=None,
        package_version=None,
        container_id="dynverse/ti_monocle_ddrtree:v0.9.9.02",
    )(
        reduction_method=reduction_method,
        max_components=max_components,
        norm_method=norm_method,
        auto_param_selection=auto_param_selection,
        filter_features=filter_features,
        filter_features_mean_expression=filter_features_mean_expression
    )
