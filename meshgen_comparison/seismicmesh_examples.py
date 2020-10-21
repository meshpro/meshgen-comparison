import numpy
import SeismicMesh

packages = [("SeismicMesh", SeismicMesh.__version__)]
colors = ("#9467bd", "#c5b0d5")  # cat20 purple


def disk(h):
    bbox = (-1.0, 1.0, -1.0, 1.0)
    domain = SeismicMesh.geometry.Disk([0.0, 0.0], 1.0)

    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox, h0=h, domain=domain, edge_length=h, verbose=0
    )
    return points, cells


def l_shape(h):
    bbox = (0.0, 1.0, 0.0, 1.0)

    bbox0 = (0.0, 1.0, 0.0, 0.5)
    rect0 = SeismicMesh.Rectangle(bbox0)

    bbox1 = (0.0, 0.5, 0.0, 1.0)
    rect1 = SeismicMesh.Rectangle(bbox1)

    corners = SeismicMesh.geometry.corners

    def union(x):
        return numpy.minimum.reduce([rect0.eval(x), rect1.eval(x)])

    pfix = numpy.vstack((corners(bbox0), corners(bbox1)))

    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox,
        domain=union,
        edge_length=h,
        pfix=pfix,
        verbose=0,
    )
    return points, cells


def rect_with_refinement(h):
    bbox = (-1.0, 1.0, -1.0, 1.0)
    domain = SeismicMesh.geometry.Rectangle([-1.0, +1.0, -1.0, +1.0])

    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox,
        h0=h,
        domain=domain,
        edge_length=lambda x: h + 0.1 * numpy.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2),
        verbose=0,
    )
    return points, cells


def ball(h):
    bbox = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)

    def ball(x):
        return numpy.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2 + x[:, 2] ** 2) - 1.0

    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox, h0=h, domain=ball, edge_length=h, verbose=0
    )
    print("b")
    points, cells = SeismicMesh.sliver_removal(
        points=points, bbox=bbox, h0=h, domain=ball, edge_length=h, verbose=False
    )
    print("c")
    return points, cells


def l_shape_3d(h):
    bbox = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

    tol = h / 10

    bbox0 = (0.0, 1.0, 0.0, 1.0, 0.0, 0.50 + tol)
    cube0 = SeismicMesh.Cube(bbox0)

    bbox1 = (0.0, 1.0, 0.5, 1.0, 0.50 - tol, 1.0)
    cube1 = SeismicMesh.Cube(bbox1)

    bbox2 = (0.5 - tol, 1.0, 0.0, 1.0, 0.50 - tol, 1.0)
    cube2 = SeismicMesh.Cube(bbox2)

    corners = SeismicMesh.geometry.corners

    def union(x):
        return numpy.minimum.reduce([cube0.eval(x), cube1.eval(x), cube2.eval(x)])

    pfix = numpy.vstack((corners(bbox0), corners(bbox1), corners(bbox2)))

    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox,
        domain=union,
        edge_length=h,
        pfix=pfix,
        verbose=0,
    )

    points, cells = SeismicMesh.sliver_removal(
        points=points,
        bbox=bbox,
        domain=union,
        edge_length=h,
        verbose=0,
    )
    return points, cells


def box_with_refinement(h):
    bbox = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)
    cube = SeismicMesh.geometry.Cube([-1.0, 1.0, -1.0, 1.0, -1.0, 1.0])

    def edge_length(x):
        return h + 0.1 * numpy.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2 + x[:, 2] ** 2)

    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox, h0=h, domain=cube, edge_length=edge_length, verbose=False
    )
    points, cells = SeismicMesh.sliver_removal(
        points=points,
        bbox=bbox,
        h0=h,
        domain=cube,
        edge_length=edge_length,
        verbose=False,
    )
    return points, cells


if __name__ == "__main__":
    import meshio

    points, cells = l_shape_3d(0.1)
    meshio.Mesh(points, {"tetra": cells}).write("out.vtk")
