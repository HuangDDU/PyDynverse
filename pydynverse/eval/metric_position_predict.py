import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from pydynverse.util.expand_matrix import expand_matrix


def calculate_position_predict(
    dataset, 
    prediction, 
    metrics=["rf_mse", "rf_rsq", "rf_nmse", "lm_mse", "lm_rsq", "lm_nmse"]
):
    """
    Compute metrics related to the prediction of the positions.
    
    Parameters:
        dataset : dict
            A dataset containing a trajectory. 必须包含至少以下键：
              - "cell_ids": 细胞标识列表。
              - "milestone_percentages": 长格式 DataFrame，包含 "cell_id", "milestone_id", "percentage"。
        prediction : dict or None
            一个预测的轨迹，其格式与 dataset 类似。
        metrics : list of str, optional
            要计算的指标，可以是 "rf_mse", "rf_rsq", "rf_nmse", "lm_mse", "lm_rsq", "lm_nmse" 中的一个或多个。
    
    Returns:
        dict: 返回一个字典，其中 summary 键下包含各指标的整体值，同时可能还包括各里程碑的详细指标。
    """
    cell_ids = dataset["cell_ids"]
    output = {"summary": {}}

    # 构造 gold_milenet_m：将 dataset["milestone_percentages"] 转换为宽格式矩阵
    gold_milenet_m = pd.pivot_table(
        dataset["milestone_percentages"],
        index="cell_id",
        columns="milestone_id",
        values="percentage",
        fill_value=0
    )
    # 使用 expand_matrix 补全（确保行顺序为 dataset["cell_ids"]）
    gold_milenet_m = expand_matrix(gold_milenet_m, rownames=cell_ids)

    # 计算 baseline_mse: 对每个里程碑（每一列），计算均方差（每个值减去该列均值的平方均值），再取所有里程碑均值
    baseline_mse = np.mean([np.mean((gold_milenet_m[col] - gold_milenet_m[col].mean())**2)
                             for col in gold_milenet_m.columns])

    if (prediction is not None and 
        len(pd.unique(prediction["milestone_percentages"]["cell_id"])) >= 3):
        # 构造预测的里程碑矩阵
        pred_milenet_m = pd.pivot_table(
            prediction["milestone_percentages"],
            index="cell_id",
            columns="milestone_id",
            values="percentage",
            fill_value=0
        )
        pred_milenet_m = expand_matrix(pred_milenet_m, rownames=cell_ids)
        # 仅保留标准差大于0的里程碑列
        cols = [col for col in pred_milenet_m.columns if pred_milenet_m[col].std() > 0]
        pred_milenet_m = pred_milenet_m[cols]

        # 随机森林模型
        if any(metric in metrics for metric in ["rf_mse", "rf_rsq", "rf_nmse"]):
            rf_mses = {}
            rf_rsqs = {}
            for col in gold_milenet_m.columns:
                # 构造数据：目标为 gold_milenet_m 某列（重命名为 "PREDICT"），预测变量为所有 pred_milenet_m 列
                target = gold_milenet_m[[col]].rename(columns={col: "PREDICT"})
                data = pd.concat([target, pred_milenet_m], axis=1)
                # 随机森林回归器：使用 5000 棵树，1 线程
                rf = RandomForestRegressor(n_estimators=5000, n_jobs=1, random_state=42)
                rf.fit(data.drop("PREDICT", axis=1), data["PREDICT"])
                preds = rf.predict(data.drop("PREDICT", axis=1))
                mse = mean_squared_error(data["PREDICT"], preds)
                rf_mses[col] = mse
                rsq = rf.score(data.drop("PREDICT", axis=1), data["PREDICT"])
                if np.isnan(rsq):
                    rsq = 1
                rf_rsqs[col] = rsq
            output["rf_mses"] = rf_mses
            output["summary"]["rf_mse"] = np.mean(list(rf_mses.values()))
            output["rf_rsqs"] = rf_rsqs
            output["summary"]["rf_rsq"] = max(np.mean(list(rf_rsqs.values())), 0)
            output["summary"]["rf_nmse"] = max(1 - output["summary"]["rf_mse"] / baseline_mse, 0)

        # 线性模型
        if any(metric in metrics for metric in ["lm_mse", "lm_rsq", "lm_nmse"]):
            lm_mses = []
            lm_rsqs = {}
            for col in gold_milenet_m.columns:
                target = gold_milenet_m[[col]].rename(columns={col: "PREDICT"})
                data = pd.concat([target, pred_milenet_m], axis=1)
                lr = LinearRegression()
                lr.fit(data.drop("PREDICT", axis=1), data["PREDICT"])
                preds = lr.predict(data.drop("PREDICT", axis=1))
                mse = np.mean((data["PREDICT"] - preds)**2)
                lm_mses.append(mse)
                rsq = lr.score(data.drop("PREDICT", axis=1), data["PREDICT"])
                lm_rsqs[col] = rsq if not np.isnan(rsq) else 1
            output["summary"]["lm_mse"] = np.mean(lm_mses)
            output["lm_rsqs"] = lm_rsqs
            output["summary"]["lm_rsq"] = max(np.mean(list(lm_rsqs.values())), 0)
            output["summary"]["lm_nmse"] = max(1 - output["summary"]["lm_mse"] / baseline_mse, 0)
    else:
        output["summary"] = {
            "rf_mse": baseline_mse,
            "rf_nmse": 0,
            "rf_rsq": 0,
            "lm_mse": baseline_mse,
            "lm_rsq": 0,
            "lm_nmse": 0
        }
    return output
