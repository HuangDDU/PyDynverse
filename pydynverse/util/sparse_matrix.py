from scipy.sparse import issparse

def is_sparse(x):
    """
    Check if an object is a sparse matrix.
    
    Parameters:
      x: Any
         对象，用于检测是否为稀疏矩阵。
    
    Returns:
      bool: 如果 x 是一个稀疏矩阵则返回 True，否则返回 False。
    """
    return issparse(x)