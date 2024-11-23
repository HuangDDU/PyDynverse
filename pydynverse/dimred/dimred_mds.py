from sklearn.metrics.pairwise import pairwise_distances
from sklearn.manifold import MDS

from .process_dimred import _process_dimred

def dimred_mds(x, ndim=2, distance_method="euclidean"):
    # TODO: 这里distance_method应该可以选择
    dis = pairwise_distances(x, metric=distance_method)
    space = MDS(n_components=ndim).fit_transform(dis)
    return _process_dimred(space)
