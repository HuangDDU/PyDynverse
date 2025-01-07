from sklearn.metrics.pairwise import pairwise_distances
from sklearn.manifold import TSNE

from .process_dimred import _process_dimred


def dimred_tsne(x, ndim=2):
    space = TSNE(n_components=ndim, init="random").fit_transform(x)
    return _process_dimred(space)
