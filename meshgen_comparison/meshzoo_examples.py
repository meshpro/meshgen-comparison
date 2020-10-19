import meshzoo

packages = [("meshzoo", meshzoo.__version__)]
colors = ("#8c564b", "#c49c94")  # cat20 brown


def disk(h):
    # tighten a bit
    h /= 1.1
    return meshzoo.disk(6, int(1 / h))
