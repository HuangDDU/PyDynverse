import matplotlib.colors as mcolors
import seaborn as sns


# RGB转化为16进制
def rgb2hex(palette):
    return [mcolors.to_hex(color) for color in palette]

# 调色板
def milestone_palette(name, n):
    if name == "auto":
        # 默认使用Set调色板
        palette = sns.color_palette("Set3")[:n]
    else:
        palette = sns.color_palette(name)[:n]

    return rgb2hex(palette)