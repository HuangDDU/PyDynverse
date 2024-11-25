from .method_choose_backend import method_choose_backend


def ti_slingshot(
    cluster_method="pam",
        ndim=20,
        shrink=1,
        reweight=True,
        reassign=True,
        thresh=0.001,
        maxit=10,
        stretch=2,
        smoother="smooth.spline",
        # shrink.method="cosine"
):
    # 选择后端执行, 备选的Docker容器ID, 传参执行
    return method_choose_backend(
        package_repository=None,
        package_name=None,
        function_name=None,
        package_version=None,
        container_id="dynverse/ti_slingshot:v1.0.3",
    )(
        cluster_method = cluster_method,
        ndim = ndim,
        shrink = shrink,
        reweight = reweight,
        reassign = reassign,
        thresh = thresh,
        maxit = maxit,
        stretch = stretch,
        smoother = smoother,
        # shrink.method = shrink.method
    )
