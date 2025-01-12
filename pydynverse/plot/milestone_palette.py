import matplotlib.colors as mcolors
import seaborn as sns

from .._settings import settings
from .._logging import logger

# RGB转化为16进制


def rgb2hex(palette):
    return [mcolors.to_hex(color) for color in palette]


def milestone_palette(name, n):
    if name == "auto":
        # 选择调色版, 颜色数量超过的话,需要自动选择
        name = settings.sns_palette
        palette = sns.color_palette(name)
        if n <= len(palette):
            palette = palette[:n]
        else:
            logger.warning(f"The number of colors({n}) is greater than the number of colors in the '{name}' palette({len(palette)}), and the 'husl' palette selection is used.")
            palette = sns.color_palette("husl", n_colors=n)
    else:
        palette = sns.color_palette(name)[:n]

    return rgb2hex(palette)
