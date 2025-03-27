import numpy as np
import pandas as pd
import scipy.sparse as sp
from pydynverse.wrap.wrap_add_trajectory import is_wrapper_with_trajectory
from pydynverse.feature.calculate_feature_importances import calculate_feature_importances
from pydynverse.feature.fi_methods import fi_ranger_rf_lite
from pydynverse.wrap.wrap_add_expression import get_expression
from pydynverse.util.expand_matrix import expand_matrix

def calculate_milestone_feature_importance(
        trajectory, 
        expression_source="expression", 
        milestones_oi=None, 
        fi_method=None, 
        verbose=False
):
#请注意，虽然这里fi_method未选取，但是其实在调用这个参数的更底层函数里指定好了默认使用的方法
    """
    目标：​计算轨迹数据中每个里程碑（Milestone）的特征重要性，即哪些基因（特征）对预测细胞（Cell）属于某个里程碑的贡献最大。其实现步骤如下：
    1.​提取表达矩阵：从轨迹对象中获取基因表达数据。
    2.​验证轨迹完整性：确保轨迹包含必要的细胞和里程碑信息。
    ​3.构建里程碑百分比矩阵：将长格式的里程碑分配数据转换为宽格式矩阵。
    ​4.计算特征重要性：调用底层函数，分析表达数据与里程碑的关系。
    ​5.重命名结果列：将输出中的 predictor_id 改为 milestone_id



    计算里程碑特征重要性得分。
    参数：
        trajectory : dict
            轨迹数据（包装器），至少包含以下关键字：
                - cell_ids"：单元格标识符列表。
                - milestone_percentages：包含cell_id、milestone_id 和 percentage。
                - milestone_id"：里程碑标识符列表。
                - 此外，expression（表达矩阵）数据可通过 get_expression 访问。
        expression_source：字符串或矩阵，可选
            用于从轨迹中提取表达式矩阵的关键字（如 “expression”）、
            或直接是表达式矩阵。
        milestones_oi : list-like, 可选
            感兴趣的里程碑。如果无，则使用轨迹[“milestone_ids”]中的所有里程碑。
        fi_method : dict, 可选
            特征重要性方法。它应该是一个字典，关键字 “fun ”指的是一个
            可调用函数，其签名为 fun(X,y,verbose)。如果为空，则使用默认方法（例如 fi_ranger_rf_lite()
            方法。
        verbose : bool, 可选
            是否打印额外信息。

    返回值
        pd.DataFrame： 一个包含三列的 DataFrame：milestone_id、feature_id 和 importance。
                      (注意：原来的 predictor_id 重命名为 milestone_id）。
    
    Calculate milestone feature importance scores.

    Parameters:
        trajectory : dict
            Trajectory data (wrapper) containing at least the following keys:
                - "cell_ids": list of cell identifiers.
                - "milestone_percentages": DataFrame in long format with columns
                  "cell_id", "milestone_id", and "percentage".
                - "milestone_ids": list of milestone identifiers.
                - Additionally, expression data is expected to be accessible via get_expression.
        expression_source : str or matrix, optional
            A key (e.g. "expression") used to extract the expression matrix from trajectory,
            or directly an expression matrix.
        milestones_oi : list-like, optional
            Milestones of interest. If None, all milestones from trajectory["milestone_ids"] are used.
        fi_method : dict, optional
            Feature importance method. It should be a dictionary with key 'fun' referring to a
            callable with signature fun(X, y, verbose). If None, a default (e.g. fi_ranger_rf_lite())
            method should be used.
        verbose : bool, optional
            Whether to print out extra information.

    Returns:
        pd.DataFrame: A DataFrame with three columns: milestone_id, feature_id, and importance.
                      (Note: the original predictor_id is renamed to milestone_id.)
    """
    
    # 提取表达矩阵，调用定义的 get_expression 接口
    expression = get_expression(trajectory, expression_source)

    # 检查 trajectory 是否包含轨迹信息（使用 is_wrapper_with_trajectory 接口）
    if not is_wrapper_with_trajectory(trajectory):
        raise ValueError("The provided trajectory does not contain trajectory information.")

    milestone_percentages = trajectory.get("milestone_percentages")
    cell_ids = trajectory.get("cell_ids")

    # 检查所有 cell_ids 是否都出现在表达矩阵的行索引中
    if not set(cell_ids).issubset(set(expression.index)):
        raise ValueError("Not all cell_ids in trajectory are present in the expression data.")
    if len(cell_ids) < 3:
        raise ValueError("Need 3 or more cells in a trajectory to determine important features.")

    # 处理里程碑：如果未指定 milestones_oi，则使用 trajectory 中的所有 milestone_ids
    if milestones_oi is None:
        milestones_oi = trajectory.get("milestone_ids")

    # 构建里程碑百分比矩阵：
    # 1. 过滤出 milestone_percentages 中 milestone_id 属于 milestones_oi 的行
    filtered = milestone_percentages[milestone_percentages["milestone_id"].isin(milestones_oi)]
    # 2. 将 long 格式转换为 wide 格式，构造行：cell_id，列：milestone_id，值：percentage，缺失值填 0
    pivot = filtered.pivot_table(index="cell_id", columns="milestone_id", 
                                 values="percentage", fill_value=0)
    # 3. 使用 expand_matrix 补全（如果原来某些 cell 没有数据，则补 0）
    milenet_m = expand_matrix(pivot, rownames=cell_ids)

    # 计算特征重要性：调用之前实现的 calculate_feature_importances
    result = calculate_feature_importances(X=expression, Y=milenet_m, fi_method=fi_method, verbose=verbose)
    # 将结果中的 predictor_id 重命名为 milestone_id
    result = result.rename(columns={"predictor_id": "milestone_id"})
    return result
