import dmsh


def disk(h):
    geo = dmsh.Circle([0.0, 0.0], 1.0)
    return dmsh.generate(geo, h)
