import pytest
import pydynverse as pdv

import matplotlib.pyplot as plt
from .test_plot_dimred import get_test_wrap_data


def test_plot_topology():
    trajectory, grouping, space = get_test_wrap_data()
    pdv.plot.plot_topology(trajectory)
    plt.savefig("/home/huang/PyCode/scRNA/PyDynverse/PyDynverse/tests/plot/img/test_plot_topology.png")


if __name__ == "__main__":
    pytest.main(["-v", __file__])
