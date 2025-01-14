import networkx as nx
from .add_milestone_coloring import add_milestone_coloring


def plot_topology(trajectory, ax=None, nx_draw_kwrags={}, layout="dot"):
    milestone_color = add_milestone_coloring(trajectory["milestone_ids"])

    milestone_network = trajectory["milestone_network"].copy()
    milestone_network["len"] = milestone_network["length"]  # 边权值对于权值的影响
    is_directed = milestone_network["directed"].any()
    G = nx.from_pandas_edgelist(milestone_network, source="from", target="to", edge_attr=["len"], create_using=nx.DiGraph if is_directed else nx.Graph)
    if layout in ["dot", "neato"]:
        pos = nx.nx_agraph.graphviz_layout(G, prog=layout)  # dot布局层级结构明显, neato布局考虑边权值
    else:
        pos = nx.spring_layout(G)

    nx.draw(G, pos, with_labels=True, node_color=[milestone_color.loc[node, "color"] for node in G.nodes], ax=ax, **nx_draw_kwrags)
