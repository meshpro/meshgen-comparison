import dmsh
import numpy

packages = [("dmsh", dmsh.__version__)]
colors = ("#2ca02c", "#98df8a")  # cat20 green


def disk(h):
    geo = dmsh.Circle([0.0, 0.0], 1.0)
    points, cells = dmsh.generate(geo, edge_size=h)
    return points, cells

    # points, cells = optimesh.cvt.quasi_newton_uniform_full(
    #     points, cells, 1.0e-2, 100, verbose=False
    # )
    # return points, cells


def rect_with_refinement(h):
    return dmsh.generate(
        dmsh.Rectangle(-1.0, 1.0, -1.0, 1.0),
        edge_size=lambda x: h + 0.1 * numpy.sqrt(x[0] ** 2 + x[1] ** 2),
        tol=1.0e-10,
    )


def l_shape(h):
    return dmsh.generate(
        dmsh.Polygon(
            [
                [-1.0, -1.0],
                [+1.0, -1.0],
                [+1.0, +0.0],
                [+0.0, +0.0],
                [+0.0, +1.0],
                [-1.0, +1.0],
            ]
        ),
        edge_size=h,
        tol=1.0e-10,
    )


if __name__ == "__main__":
    import meshio

    points, cells = disk(0.1)
    meshio.Mesh(points, {"triangle": cells}).write("out.vtk")
