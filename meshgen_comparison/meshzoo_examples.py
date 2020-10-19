import meshzoo

packages = [("meshzoo", meshzoo.__version__)]
colors = ("#2ca02c", "#98df8a")


def disk(h):
    # tighten a bit
    h /= 1.1
    return meshzoo.disk(6, int(1 / h))
