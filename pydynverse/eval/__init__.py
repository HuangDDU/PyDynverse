import os.path
import pandas as pd
from .calculate_metrics import calculate_metrics
from .metric_isomorphic import calc_isomorphic


metrics = pd.read_csv(f"{os.path.dirname(__file__)}/metrics.csv", sep='\t')

__all__ = [
    "metrics",
    "calculate_metrics",
    "calc_isomorphic"
]
