import matplotlib.colors as mcolors


def mix_colors(milid, milpct, milestone_colors):
    # 按照percentage对于混合颜色后，转化为16进制
    # milid: milestone_id, milpct: milestone_percentange
    color_rgb = milestone_colors.loc[milid].apply(lambda x: (x.array*milpct.array).sum())  # 分别对RGB三个通道按照percentage加权
    # 暂时不用阈值设置，加权之后肯定在0-1之间
    # color_rgb[color_rgb < 0] = 0
    # color_rgb[color_rgb > 1] = 1
    return mcolors.to_hex(color_rgb.array)
