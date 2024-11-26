from .method_choose_backend import method_choose_backend


def ti_grandprix(
    n_inducing_points=4,
    latent_prior_var=0.1,
    latent_var=0.028
):
    return method_choose_backend(
        package_repository=None,
        package_name=None,
        function_name=None,
        package_version=None,
        container_id="dynverse/ti_grandprix:v0.9.9.01"
    )(
        n_inducing_points=n_inducing_points,
        latent_prior_var=latent_prior_var,
        latent_var=latent_var
    )
