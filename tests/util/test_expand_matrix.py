import numpy as np
import pandas as pd
import pytest
from pydynverse.util.expand_matrix import expand_matrix  

def test_expand_matrix():
    # 构造一个 3x4 的 DataFrame，行标签为 ["a", "c", "d"]，列标签为 ["D", "F", "H", "I"]
    data = np.arange(12).reshape(3, 4)
    df = pd.DataFrame(data, index=["a", "c", "d"], columns=["D", "F", "H", "I"])
    
    # 期望的行标签和列标签
    desired_rows = ["a", "b", "c", "d", "e"]  # 其中 "b" 和 "e" 原始中没有
    desired_cols = list("ABCDEFGHIJ")[:10]      # 例如前 10 个大写字母
    
    # 调用 expand_matrix，缺失位置用 0 填充
    expanded = expand_matrix(df, rownames=desired_rows, colnames=desired_cols, fill=0)
    
    # 检查输出矩阵尺寸
    assert expanded.shape == (len(desired_rows), len(desired_cols))
    
    # 对于原矩阵中存在的行和列，值应保持不变
    for row in ["a", "c", "d"]:
        for col in df.columns:
            if col in desired_cols:
                assert expanded.loc[row, col] == df.loc[row, col]
    
    # 对于原矩阵中不存在的行（如 "b", "e"）其所有值应为填充值 0
    for row in ["b", "e"]:
        for col in desired_cols:
            assert expanded.loc[row, col] == 0

if __name__ == "__main__":
    pytest.main(["-v", __file__])