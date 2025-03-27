import numpy as np
import pandas as pd
import pytest
from scipy.stats import ks_2samp, ranksums

from pydynverse.feature.fi_methods import fi_ranger_rf_lite  # 默认方法，但用 dummy_fi_method 替代
from pydynverse.feature.calculate_overall_feature_importance import calculate_overall_feature_importance
from pydynverse.eval.metric_featureimp import (
    calculate_featureimp_cor,
    calculate_featureimp_enrichment,
    _calculate_featureimp_cor  # 内部函数测试
)


def dummy_fi_method():
    """
    一个简单的 dummy 特征重要性方法：
    对于给定的 X, y，直接返回 X 各列的样本方差（使用 ddof=1），结果以字典形式返回。
    """
    def fun(X, y, verbose=False):
        return X.var(ddof=1).to_dict()
    return {'fun': fun}


# 构造测试轨迹（trajectory）数据

@pytest.fixture
def dummy_dataset():
    """
    构造一个 dataset 轨迹，包含：
      - cell_ids: 四个细胞
      - milestone_percentages: 每个细胞对应一个里程碑（前两细胞为 m1，后两为 m2）
      - milestone_ids: ["m1", "m2"]
      - expression: 一个 4x5 的表达矩阵（行为细胞，列为基因），取值固定以便计算方差
      - prior_information: 包含先验特征列表，用于 enrichment 测试
    """
    cell_ids = ["c1", "c2", "c3", "c4"]
    # 构造表达矩阵，每个基因在 4 个细胞中的取值确定（便于计算方差）
    # gene0: [1,2,3,4], gene1: [2,2,2,2], gene2: [1,3,5,7], gene3: [4,5,6,7], gene4: [1,1,1,1]
    expr = pd.DataFrame({
        "gene0": [1, 2, 3, 4],
        "gene1": [2, 2, 2, 2],
        "gene2": [1, 3, 5, 7],
        "gene3": [4, 5, 6, 7],
        "gene4": [1, 1, 1, 1]
    }, index=cell_ids)
    # milestone_percentages：简单构造，每个细胞只有一个里程碑分配（以百分比 1.0 表示）
    milestone_percentages = pd.DataFrame({
        "cell_id": cell_ids,
        "milestone_id": ["m1", "m1", "m2", "m2"],
        "percentage": [1.0, 1.0, 1.0, 1.0]
    })
    return {
        "cell_ids": cell_ids,
        "milestone_percentages": milestone_percentages,
        "milestone_ids": ["m1", "m2"],
        "pydynwrap:with_trajectory": True,
        "expression": expr,
        # 先验信息，用于 enrichment：假设 dataset 中认为重要的特征为 gene0 和 gene2
        "prior_information": {"features_id": ["gene0", "gene2"]}
    }

@pytest.fixture
def dummy_prediction():
    """
    构造一个 prediction 轨迹，与 dataset 类似，但表达矩阵做一个线性变换（例如乘以2）
    使得两个轨迹的整体特征重要性存在严格的线性关系。
    """
    cell_ids = ["c1", "c2", "c3", "c4"]
    # 使用 dataset 表达矩阵的2倍
    expr = pd.DataFrame({
        "gene0": [2, 4, 6, 8],
        "gene1": [4, 4, 4, 4],
        "gene2": [2, 6, 10, 14],
        "gene3": [8, 10, 12, 14],
        "gene4": [2, 2, 2, 2]
    }, index=cell_ids)
    milestone_percentages = pd.DataFrame({
        "cell_id": cell_ids,
        "milestone_id": ["m1", "m1", "m2", "m2"],
        "percentage": [1.0, 1.0, 1.0, 1.0]
    })
    return {
        "cell_ids": cell_ids,
        "milestone_percentages": milestone_percentages,
        "milestone_ids": ["m1", "m2"],
        "pydynwrap:with_trajectory": True,
        "expression": expr
    }


# 计算整体特征重要性时，使用 dummy_fi_method 以获得可预期结果
# 对于 dataset 表达矩阵，按列计算样本方差（ddof=1）：
# gene0: variance([1,2,3,4]) = 1.6667, gene1: variance([2,2,2,2]) = 0, gene2: variance([1,3,5,7]) = 6.6667,
# gene3: variance([4,5,6,7]) = 1.6667, gene4: variance([1,1,1,1]) = 0.
# 对于 prediction 表达矩阵（是 dataset 的2倍），方差会放大 4 倍：
# gene0: 1.6667*4 = 6.6667, gene1: 0, gene2: 6.6667*4 = 26.6667, gene3: 1.6667*4 = 6.6667, gene4: 0.
# calculate_overall_feature_importance 仅做 groupby 平均（这里每个细胞只有一种预测，所以结果即为上面的值）。
# 测试 _calculate_featureimp_cor（内部函数）

def test__calculate_featureimp_cor_normal():
    # 构造两组整体特征重要性 DataFrame（来自 dataset 和 prediction）
    dataset_imp = pd.DataFrame({
        "feature_id": ["gene0", "gene1", "gene2", "gene3", "gene4"],
        "importance": [1.6667, 0, 6.6667, 1.6667, 0]
    })
    pred_imp = pd.DataFrame({
        "feature_id": ["gene0", "gene1", "gene2", "gene3", "gene4"],
        "importance": [6.6667, 0, 26.6667, 6.6667, 0]
    })
    result = _calculate_featureimp_cor(dataset_imp, pred_imp)
    # 由于两者严格存在线性关系（gene1和gene4均为0），
    # 计算时整体标准差不为0（因为 gene0, gene2, gene3 有变化），相关性应为1。
    assert pytest.approx(result["featureimp_cor"], rel=1e-3) == 1
    assert pytest.approx(result["featureimp_wcor"], rel=1e-3) == 1

def test__calculate_featureimp_cor_constant():
    # 如果 dataset_imp 所有值都相同，则返回0
    dataset_imp = pd.DataFrame({
        "feature_id": ["gene0", "gene1", "gene2"],
        "importance": [0.5, 0.5, 0.5]
    })
    pred_imp = pd.DataFrame({
        "feature_id": ["gene0", "gene1", "gene2"],
        "importance": [6.6667, 0, 26.6667]
    })
    result = _calculate_featureimp_cor(dataset_imp, pred_imp)
    assert result["featureimp_cor"] == 0
    assert result["featureimp_wcor"] == 0

def test__calculate_featureimp_cor_missing_features():
    # 部分特征缺失，通过 outer join 补0后仍能计算相关性
    dataset_imp = pd.DataFrame({
        "feature_id": ["gene0", "gene1"],
        "importance": [1.6667, 0]
    })
    pred_imp = pd.DataFrame({
        "feature_id": ["gene1", "gene2"],
        "importance": [0, 26.6667]
    })
    result = _calculate_featureimp_cor(dataset_imp, pred_imp)
    # 此时 join 后 gene0: (1.6667,0), gene1: (0,0), gene2: (0,26.6667)
    # 只要标准差不为0（整体数据有变化），返回结果应不为负
    assert result["featureimp_cor"] >= 0
    assert result["featureimp_wcor"] >= 0


# 测试 calculate_featureimp_cor 整体函数（依赖 calculate_overall_feature_importance 实际调用）

def test_calculate_featureimp_cor_integration(dummy_dataset, dummy_prediction):
    # 使用 dummy_fi_method 替代默认 fi_method，以便得到可预测的整体重要性结果
    result = calculate_featureimp_cor(
        dataset=dummy_dataset,
        prediction=dummy_prediction,
        expression_source="expression",
        fi_method=dummy_fi_method()
    )
    # 根据上面说明：
    # dummy_dataset overall importance：gene0=1.6667, gene1=0, gene2=6.6667, gene3=1.6667, gene4=0
    # dummy_prediction overall importance：gene0=6.6667, gene1=0, gene2=26.6667, gene3=6.6667, gene4=0
    # 经过 _calculate_featureimp_cor，相关性应为 1。
    assert pytest.approx(result["featureimp_cor"], rel=1e-3) == 1
    assert pytest.approx(result["featureimp_wcor"], rel=1e-3) == 1

def test_calculate_featureimp_cor_invalid(dummy_dataset):
    # 如果 prediction 中 milestone_percentages 的 cell_id 数量不足（小于3），直接返回 0
    invalid_prediction = {
        "milestone_percentages": pd.DataFrame({
            "cell_id": ["c1", "c1"],  # 重复值，唯一值数量不足3
            "milestone_id": ["m1", "m1"],
            "percentage": [1.0, 1.0]
        })
    }
    result = calculate_featureimp_cor(dummy_dataset, invalid_prediction)
    assert result == {'featureimp_cor': 0, 'featureimp_wcor': 0}


# 测试 calculate_featureimp_enrichment 整体函数

def test_calculate_featureimp_enrichment_not_enough_notsel(dummy_dataset, dummy_prediction):
    """
    情形1：若 prediction 得到的整体重要性结果中，未选中（not sel）的特征数不超过2，
    则应返回 {featureimp_ks: 1, featureimp_wilcox: 1}.
    为此我们构造一个预测轨迹，其表达矩阵只产生2个基因的重要性结果。
    """
    # 修改 dummy_prediction 的 expression，使得只有 gene0 和 gene1 出现
    cell_ids = dummy_prediction["cell_ids"] if "cell_ids" in dummy_prediction else ["c1", "c2", "c3", "c4"]
    expr = pd.DataFrame({
        "gene0": [2,4,6,8],
        "gene1": [4,4,4,4]
    }, index=cell_ids)
    pred = dummy_prediction.copy()
    pred["expression"] = expr
    # 此时 calculate_overall_feature_importance 计算得到的整体重要性只有 gene0 和 gene1
    result = calculate_featureimp_enrichment(
        dataset=dummy_dataset,
        prediction=pred,
        expression_source="expression",
        fi_method=dummy_fi_method()
    )
    assert result == {'featureimp_ks': 1, 'featureimp_wilcox': 1}

def test_calculate_featureimp_enrichment_valid(dummy_dataset, dummy_prediction):
    """
    情形2：若未选中特征数足够 (>2)，则应进行 KS 和 Wilcoxon 检验，
    返回的 p 值在 [0,1] 范围内。由于 dummy_fi_method 返回的整体重要性结果相同（两轨迹均相同），
    KS 检验的 p-value 应为 1，Wilcoxon 检验 p-value 也为 1，从而返回 {ks:1, wilcox: 0}.
    """
    # 修改 dataset 的先验信息，使其仅包含部分基因
    dummy_dataset["prior_information"]["features_id"] = ["gene0", "gene2"]
    result = calculate_featureimp_enrichment(
        dataset=dummy_dataset,
        prediction=dummy_prediction,
        expression_source="expression",
        fi_method=dummy_fi_method()
    )
    # 计算 overall importance 对于 dummy_dataset 与 dummy_prediction（已在前面计算）
    # 结果为：每个基因 importance：gene0=1.6667, gene1=0, gene2=6.6667, gene3=1.6667, gene4=0
    # sel 为 [gene0, gene2] => [1.6667, 6.6667]； notsel 为 [gene1, gene3, gene4] => [0, 1.6667, 0]
    # 对于完全相同的分布， KS 检验 p-value 应为 1，ranksums p-value 应为 1，
    # 则返回 {featureimp_ks: 1, featureimp_wilcox: 1-1=0}.
    assert 0 <= result['featureimp_ks'] <= 1
    assert 0 <= result['featureimp_wilcox'] <= 1
    # 由于分布极为相似，期望 ks 接近1，wilcox 接近0
    assert pytest.approx(result['featureimp_ks'], rel=1e-3) == 1
    assert pytest.approx(result['featureimp_wilcox'], rel=1e-3) == 0

def test_calculate_featureimp_enrichment_invalid(dummy_dataset):
    # 如果 prediction 无效（cell_id 数量不足），应返回 {0,0}
    invalid_prediction = {
        "milestone_percentages": pd.DataFrame({
            "cell_id": ["c1", "c1"],
            "milestone_id": ["m1", "m1"],
            "percentage": [1.0, 1.0]
        })
    }
    result = calculate_featureimp_enrichment(dummy_dataset, invalid_prediction)
    assert result == {'featureimp_ks': 0, 'featureimp_wilcox': 0}

if __name__ == "__main__":
    pytest.main(["-v", __file__])
