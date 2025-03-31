import pandas as pd


topologies_with_same_n_milestones = {
    "linear": pd.DataFrame({
        "from": ["A", "B", "C", "D", "E"],
        "to": ["B", "C", "D", "E", "F"],
    }),
    "bifurcation": pd.DataFrame({
        "from": ["A", "B", "B", "C", "D"],
        "to": ["B", "C", "D", "E", "F"],
    }),
    "multifurcating": pd.DataFrame({
        "from": ["A", "B", "B", "B", "C"],
        "to": ["B", "C", "D", "E", "F"],
    }),
    "tree": pd.DataFrame({
        "from": ["A", "B", "B", "C", "C"],
        "to": ["B", "C", "D", "E", "F"],
    }),
    "cycle": pd.DataFrame({
        "from": ["A", "B", "C", "D", "E", "F"],
        "to": ["B", "C", "D", "E", "F", "A"],
    }),
    "connected": pd.DataFrame({
        "from": ["A", "B", "B", "D", "E", "C"],
        "to": ["B", "C", "D", "E", "A", "F"],
    }),
    "disconnected": pd.DataFrame({
        "from": ["A", "B", "C", "D", "E"],
        "to": ["B", "C", "A", "E", "F"],
    })
}
