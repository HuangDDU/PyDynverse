import os
import sys


def pytest_configure(config):
    # pytest配置: 测试工作路径在项目文件夹下, 只对直接运行有效
    # 其他方式1: .bashrc文件夹下添加, 可以在任意路径下访问
    # 其他方式2: pytest.ini文件夹
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    sys.path.append(project_root)
