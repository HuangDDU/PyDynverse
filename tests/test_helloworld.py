import pytest

def test_helloworld():
    s = "hello world"
    assert s == "hello world", "hello world failed!"

if __name__ == "__main__":
    pytest.main() # 会测试当前目录下所有文件