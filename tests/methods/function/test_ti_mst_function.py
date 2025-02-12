import pytest
import pydynverse as pdv


def test_ti_mst_function():

    pdv.settings.backend = "python_function"
    dataset = pdv.data.load_simulation_data()

    method = pdv.methods.ti_mst()
    dataset = pdv.wrap.infer_trajectory(dataset, method)

    assert "progressions" in dataset


if __name__ == "__main__":
    pytest.main(["-W", "ignore"])
