import numpy
import SeismicMesh

packages = [("SeismicMesh", SeismicMesh.__version__)]
colors = ("#9467bd", "#c5b0d5")  # cat20 purple


def disk(h):
    domain = SeismicMesh.geometry.Disk([0.0, 0.0], 1.0)

    points, cells = SeismicMesh.generate_mesh(
        h0=h, domain=domain, edge_length=h, verbose=0
    )
    return points, cells


def l_shape(h):
    rect0 = SeismicMesh.Rectangle([-1.0, 1.0, -1.0, 0.0])
    rect1 = SeismicMesh.Rectangle([-1.0, 0.0, -1.0, 1.0])
    domain = SeismicMesh.Union([rect0, rect1])
    points, cells = SeismicMesh.generate_mesh(domain=domain, edge_length=h, verbose=0)
    return points, cells


def rect_with_refinement(h):
    domain = SeismicMesh.geometry.Rectangle((-1.0, +1.0, -1.0, +1.0))

    def f(x):
        return h + 0.1 * numpy.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2)

    points, cells = SeismicMesh.generate_mesh(
        domain=domain, h0=h, edge_length=f, verbose=0
    )
    return points, cells


def quarter_annulus(h):
    rect = SeismicMesh.Rectangle((0.0, 1.0, 0.0, 1.0))
    disk0 = SeismicMesh.Disk([0.0, 0.0], 0.25)
    diff0 = SeismicMesh.Difference([rect, disk0])

    disk1 = SeismicMesh.Disk([0.0, 0.0], 1.0)
    quarter = SeismicMesh.Intersection([diff0, disk1])

    points, cells = SeismicMesh.generate_mesh(
        domain=quarter,
        edge_length=lambda x: h + 0.1 * numpy.abs(disk0.eval(x)),
        h0=h,
        verbose=0,
    )
    return points, cells


def ball(h):
    ball = SeismicMesh.Ball([0.0, 0.0, 0.0], 1.0)
    points, cells = SeismicMesh.generate_mesh(
        domain=ball, h0=h, edge_length=h, verbose=0
    )
    points, cells = SeismicMesh.sliver_removal(
        domain=ball, points=points, h0=h, edge_length=h, verbose=0
    )
    return points, cells


def cylinder(h):
    cylinder = SeismicMesh.Cylinder(h=1.0, r=1.0)
    points, cells = SeismicMesh.generate_mesh(domain=cylinder, edge_length=h, verbose=0)
    points, cells = SeismicMesh.sliver_removal(
        points=points, domain=cylinder, edge_length=h, verbose=0
    )
    return points, cells


def l_shape_3d(h):
    tol = h / 10

    cube0 = SeismicMesh.Cube((-1.0, 1.0, -1.0, 1.0, -1.0, tol))
    cube1 = SeismicMesh.Cube((-1.0, 1.0, 0.0, 1.0, -tol, 1.0))
    cube2 = SeismicMesh.Cube((-tol, 1.0, -1.0, 1.0, -tol, 1.0))
    domain = SeismicMesh.Union([cube0, cube1, cube2])

    points, cells = SeismicMesh.generate_mesh(domain=domain, edge_length=h, verbose=0)
    points, cells = SeismicMesh.sliver_removal(
        domain=domain, points=points, edge_length=h, verbose=0
    )
    return points, cells


def box_with_refinement(h):
    cube = SeismicMesh.geometry.Cube((-1.0, 1.0, -1.0, 1.0, -1.0, 1.0))

    def edge_length(x):
        return h + 0.1 * numpy.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2 + x[:, 2] ** 2)

    points, cells = SeismicMesh.generate_mesh(
        domain=cube, h0=h, edge_length=edge_length, verbose=0
    )
    points, cells = SeismicMesh.sliver_removal(
        domain=cube,
        points=points,
        h0=h,
        edge_length=edge_length,
        verbose=0,
    )
    return points, cells


if __name__ == "__main__":
    import meshio

    points, cells = quarter_annulus(0.01)
    meshio.Mesh(points, {"triangle": cells}).write("out.vtk")
