import pytest
from pydynverse.util import random_time_string
import datetime
import re

def test_random_time_string_without_name():
    # 调用函数
    time_string = random_time_string()

    # 检查时间部分
    assert re.match(r'\d{8}_\d{6}__', time_string), "时间部分格式不正确"

    # 检查随机字符部分
    parts = time_string.split('__')
    assert len(parts) == 2, "字符串分割后应有两个部分"
    assert len(parts[1]) == 10, "随机字符部分长度不为10"
    assert all(c.isalnum() for c in parts[1]), "随机字符部分包含非法字符"

def test_random_time_string_with_name():
    # 定义一个名称
    name = "test_name"
    # 调用函数
    time_string = random_time_string(name)

    # 检查时间部分
    assert re.match(r'\d{8}_\d{6}__' + re.escape(name) + '__', time_string), "时间部分和名称部分格式不正确"

    # 检查随机字符部分
    parts = time_string.split('__')
    assert len(parts) == 3, "字符串分割后应有三个部分"
    assert len(parts[2]) == 10, "随机字符部分长度不为10"
    assert all(c.isalnum() for c in parts[2]), "随机字符部分包含非法字符"

if __name__ == "__main__":
    pytest.main(["-v", __file__])
