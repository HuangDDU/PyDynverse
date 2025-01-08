import pytest
import pydynverse as pdv

import matplotlib.pyplot as plt
from .test_plot_dimred import get_test_wrap_data


def test_plot_graph():
    trajectory, grouping, space = get_test_wrap_data()
    pdv.plot.plot_graph(trajectory)
    plt.savefig("/home/huang/PyCode/scRNA/PyDynverse/PyDynverse/tests/plot/img/test_plot_graph.png")


if __name__ == "__main__":
    pytest.main(["-v", __file__])
