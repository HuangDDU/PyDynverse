from .method_choose_backend import method_choose_backend


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
        package_repository=None,
        package_name=None,
        function_name=None,
        package_version=None,
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


def ti_paga_wrapper():
    # 脚本方式调用paga的封装
    pass