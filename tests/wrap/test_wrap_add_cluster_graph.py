import pytest
import pydynverse as pdv

from .test_wrap_add_waypoints import get_test_wrap_data


def test_add_cluster_graph():
    test_wrap_data = get_test_wrap_data()
    dataset = test_wrap_data["dataset"]
    milestone_network = test_wrap_data["milestone_network"]
    grouping = ["W", "X", "X", "Z", "Z", "Z"]

    trajectory = pdv.wrap.add_cluster_graph(dataset, milestone_network, grouping)

    assert trajectory["milestone_percentages"].query("`percentage`==1")["milestone_id"].tolist() == grouping




if __name__ == "__main__":
    pytest.main(["-v", __file__])
