from sklearn.manifold import Isomap

from .process_dimred import _process_dimred


def dimred_iosmap(x, ndim=2):

    space = Isomap(n_components=2).fit_transform(x)
    return _process_dimred(space)
