import datetime
import random
import string

# from .._logging import logger

def random_time_string(name=None):
    # 获取当前时间并格式化
    current_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    # 生成随机字符串
    random_chars = ''.join(random.choices(
        string.ascii_letters + string.digits, k=10))
    # 构建最终的字符串
    if name:
        time_string = f"{current_time}__{name}__{random_chars}"
    else:
        time_string =  f"{current_time}__{random_chars}"

    # logger.debug(f"Time String: {time_string}")
    return time_string


if __name__ == "__main__":
    # logger.setLevel("DEBUG")
    print(random_time_string()) # 不带名字
    print(random_time_string("example")) # 不带名字