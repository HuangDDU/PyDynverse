import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, ranksums
from pydynverse.feature.fi_methods import fi_ranger_rf_lite
from pydynverse.feature.calculate_overall_feature_importance import calculate_overall_feature_importance

def calculate_featureimp_cor(
        dataset, 
        prediction, 
        expression_source=None, #这里作者指定的是dataset["expression_source"]
        fi_method=fi_ranger_rf_lite()
    ):
    """
    Compare feature importances derived by both trajectories.
    
    Parameters:
        dataset : dict
            A dataset (trajectory) containing expression data and trajectory information.
        prediction : dict
            A predicted trajectory.
        expression_source : str, matrix or callable, optional
            The expression data matrix or key to obtain it.
        fi_method : dict, optional
            A feature importance method (default uses fi_ranger_rf_lite()).
    
    Returns:
        dict: A dictionary with keys:
            - 'featureimp_cor': the correlation between feature importances,
            - 'featureimp_wcor': a weighted correlation.
        If prediction is None or does not contain at least 3 cells in milestone_percentages,
        returns {'featureimp_cor': 0, 'featureimp_wcor': 0}.
    """
    # 检查 prediction 中的里程碑信息是否足够
    if (prediction is not None and 
        len(prediction['milestone_percentages']['cell_id'].unique()) >= 3):
        
        dataset_imp = calculate_overall_feature_importance(
            trajectory=dataset,
            expression_source=expression_source,
            fi_method=fi_method
        )
        pred_imp = calculate_overall_feature_importance(
            trajectory=prediction,
            expression_source=expression_source,
            fi_method=fi_method
        )
        return _calculate_featureimp_cor(dataset_imp, pred_imp)
    else:
        return {'featureimp_cor': 0, 'featureimp_wcor': 0}


def _calculate_featureimp_cor(dataset_imp, pred_imp):
    """
    内部函数：计算两个特征重要性数据框之间的相关性和加权相关性。
    
    Parameters:
        dataset_imp : pd.DataFrame
            DataFrame 包含两列：'feature_id' 和 'importance' (来自 dataset)
        pred_imp : pd.DataFrame
            DataFrame 包含两列：'feature_id' 和 'importance' (来自 prediction)
    
    Returns:
        dict: {'featureimp_cor': float, 'featureimp_wcor': float}
    """
    # 使用 outer join 将两个数据框合并，并用0填充缺失值
    join = pd.merge(
        dataset_imp.rename(columns={'importance': 'dataset_imp'}),
        pred_imp.rename(columns={'importance': 'pred_imp'}),
        on='feature_id',
        how='outer'
    ).fillna(0)
    
    # 如果其中一个数据集的标准差为0，则返回0
    if join['dataset_imp'].std() == 0 or join['pred_imp'].std() == 0:
        return {'featureimp_cor': 0, 'featureimp_wcor': 0}
    
    # 计算皮尔逊相关系数，并取不小于0的值
    featureimp_cor = max(join['dataset_imp'].corr(join['pred_imp']), 0)
    
    # 计算加权协方差矩阵，权重取 dataset_imp
    cov_wt = np.cov(join['dataset_imp'], join['pred_imp'], aweights=join['dataset_imp'])
    # 根据协方差矩阵计算相关系数
    featureimp_wcor = max(cov_wt[0, 1] / np.sqrt(cov_wt[0, 0] * cov_wt[1, 1]), 0)
    
    return {'featureimp_cor': featureimp_cor, 'featureimp_wcor': featureimp_wcor}


def calculate_featureimp_enrichment(
        dataset,
        prediction,
        expression_source=None, 
        fi_method=fi_ranger_rf_lite()
    ):
    """
    Compare enrichment in finding back the most important genes.
    
    Parameters:
        dataset : dict
            A dataset containing expression data, trajectory, and prior information.
        prediction : dict
            A predicted trajectory.
        expression_source : str, matrix or callable, optional
            The expression data matrix or key.
        fi_method : dict, optional
            A feature importance method (default uses fi_ranger_rf_lite()).
    
    Returns:
        dict: A dictionary with keys:
            - 'featureimp_ks': p-value from the KS test,
            - 'featureimp_wilcox': adjusted p-value (1 - p) from the Wilcoxon test.
        如果 prediction 不符合条件，则返回 {'featureimp_ks': 0, 'featureimp_wilcox': 0}.
    """
    try:
        # 判断 prediction 中是否有足够的里程碑细胞
        if (prediction is not None and 
            len(prediction['milestone_percentages']['cell_id'].unique()) >= 3):
            
            pred_imp = calculate_overall_feature_importance(
                trajectory=prediction,
                expression_source=expression_source,
                fi_method=fi_method
            )
            
            # 获取 dataset 中预先定义的重要特征
            dataset_features = dataset['prior_information']['features_id']
            
            # 分别提取属于和不属于预定义特征的 importance 值
            sel = pred_imp.loc[pred_imp['feature_id'].isin(dataset_features), 'importance']
            notsel = pred_imp.loc[~pred_imp['feature_id'].isin(dataset_features), 'importance']
            
            if len(notsel) > 2:
                ks = ks_2samp(sel, notsel, alternative='greater')
                wilcox = ranksums(sel, notsel, alternative='greater')
                return {
                    'featureimp_ks': ks.pvalue,
                    'featureimp_wilcox': 1 - wilcox.pvalue
                }
            else:
                return {'featureimp_ks': 1, 'featureimp_wilcox': 1}
        else:
            return {'featureimp_ks': 0, 'featureimp_wilcox': 0}
    except Exception as e:
        print("featureimp_enrichment errored! check reason!", e)
        return {'featureimp_ks': 0, 'featureimp_wilcox': 0}
