import meshzoo

desc = f"meshzoo {meshzoo.__version__}"


def disk(h):
    # tighten a bit
    h /= 1.1
    return meshzoo.disk(6, int(1 / h))
