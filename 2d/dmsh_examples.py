import dmsh
# import optimesh


def disk(h):
    geo = dmsh.Circle([0.0, 0.0], 1.0)
    points, cells = dmsh.generate(geo, edge_size=h)
    return points, cells

    # points, cells = optimesh.cvt.quasi_newton_uniform_full(
    #     points, cells, 1.0e-2, 100, verbose=False
    # )
    # return points, cells
