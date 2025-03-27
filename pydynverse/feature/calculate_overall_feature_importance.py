from pydynverse.feature.fi_methods import  fi_ranger_rf_lite
from pydynverse.feature.calculate_milestone_feature_importance import calculate_milestone_feature_importance
import pandas as pd

def calculate_overall_feature_importance(
        trajectory, 
        expression_source="expression", 
        fi_method=None, 
        verbose=False
        ):
    """
    Calculating overall feature importances across a trajectory.
    
    This function calls `calculate_milestone_feature_importance` to compute feature importance
    for individual milestones (e.g. branching points), then aggregates the results across milestones
    by averaging the importance values for each feature.
    
    Parameters:
        trajectory : dict
            A trajectory object containing expression data and trajectory information.
        expression_source : str, matrix or callable, default "expression"
            The expression matrix to use. If a matrix is provided, it is used directly.
            If a string is provided, then trajectory[expression_source] is used.
            If a callable is provided, it will be called to obtain the expression.
        fi_method : dict, optional
            A feature importance method. If None, the default `fi_ranger_rf_lite()` is used.
        verbose : bool, optional
            Whether to print out extra information.
    
    Returns:
        pd.DataFrame:
            A DataFrame with at least two columns: 'feature_id' and 'importance', where
            'importance' is the mean importance across milestones, sorted in descending order.
    """
    # 如果未指定特征重要性方法，则导入默认方法（fi_ranger_rf_lite）
    if fi_method is None:
        fi_method = fi_ranger_rf_lite()
    
    # 调用 milestone 级别的特征重要性计算
    milestone_importances = calculate_milestone_feature_importance(
        trajectory=trajectory,
        expression_source=expression_source,
        fi_method=fi_method,
        verbose=verbose
    )
    
    # 聚合：按 feature_id 分组并计算平均重要性，然后按降序排序
    overall_importance = (
        milestone_importances
        .groupby("feature_id", as_index=False)["importance"]
        .mean()
        .sort_values(by="importance", ascending=False)
    )
    
    return overall_importance
