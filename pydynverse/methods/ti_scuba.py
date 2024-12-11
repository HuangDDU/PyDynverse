from .method_choose_backend import method_choose_backend


def ti_scuba(
    rigorous_gap_stats=True,
    N_dim=2,
    low_gene_threshold=1,
    low_gene_fraction_max=0.7,
    min_split=15,
    min_percentage_split=0.25
):
    return method_choose_backend(
        package_repository=None,
        package_name=None,
        function_name=None,
        package_version=None,
        container_id="dynverse/ti_scuba:v0.9.9.01",
    )(
        rigorous_gap_stats=rigorous_gap_stats,
        N_dim=N_dim,
        low_gene_threshold=low_gene_threshold,
        low_gene_fraction_max=low_gene_fraction_max,
        min_split=min_split,
        min_percentage_split=min_percentage_split
    )
