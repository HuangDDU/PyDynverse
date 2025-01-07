import pytest
import pydynverse as pdv
from scipy.sparse import csc_matrix


def test_load_simulation_data():
    dataset = pdv.data.load_simulation_data()
    assert isinstance(dataset["counts"], csc_matrix), "dataset['counts'] is not a csc_matrix"
    assert isinstance(dataset["expression"], csc_matrix), "dataset['expression'] is not a csc_matrix"


if __name__ == "__main__":
    pytest.main(["-v", __file__])
