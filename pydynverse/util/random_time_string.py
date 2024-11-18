import datetime
import random
import string

def random_time_string(name=None):
    # 获取当前时间并格式化
    current_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # 生成随机字符串
    random_chars = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
    
    # 构建最终的字符串
    if name:
        return f"{current_time}__{name}__{random_chars}"
    else:
        return f"{current_time}__{random_chars}"

if __name__ == "__main__":
    # 示例调用
    print(random_time_string())  # 不带名字
    print(random_time_string("example"))  # 带名字