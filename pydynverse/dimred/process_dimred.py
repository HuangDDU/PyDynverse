import pandas as pd

def _process_dimred(space):
  space = pd.DataFrame(space)
  space.columns = [f"comp_{i}" for i in range(space.shape[1])] # 添加降维后的维度名称
  return space
