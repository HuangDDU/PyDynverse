import networkx as nx
from .add_milestone_coloring import add_milestone_coloring


def plot_topology(trajectory, ax=None, nx_draw_kwrags={}):
    milestone_color = add_milestone_coloring(trajectory["milestone_ids"])

    milestone_network = trajectory["milestone_network"]
    is_directed = milestone_network["directed"].any()
    G = nx.from_pandas_edgelist(milestone_network, source="from", target="to", edge_attr=True, create_using=nx.DiGraph if is_directed else nx.Graph)
    pos = nx.nx_agraph.graphviz_layout(G, prog="dot")  # graphviz的dot有向无环图布局更加合适

    nx.draw(G, pos, with_labels=True, node_color=[milestone_color.loc[node, "color"] for node in G.nodes], ax=ax, **nx_draw_kwrags)
