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


def rect_with_refinement(h):
    bbox = (-1.0, 1.0, -1.0, 1.0)
    domain = SeismicMesh.geometry.Rectangle([-1.0, +1.0, -1.0, +1.0])

    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox,
        h0=h,
        domain=domain,
        edge_length=lambda x: h + 0.1 * numpy.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2),
    )
    return points, cells


def ball(h):
    bbox = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)

    def ball(x):
        return numpy.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2 + x[:, 2] ** 2) - 1.0

    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox, h0=h, domain=ball, edge_length=h, verbose=False
    )
    points, cells = SeismicMesh.sliver_removal(
        points=points, bbox=bbox, h0=h, domain=ball, edge_length=h, verbose=False
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

    # points, cells = rect_with_refinement(0.1)
    points, cells = ball(0.1)
    # points, cells = box_with_refinement(0.1)
    meshio.Mesh(points, {"tetra": cells}).write("out.vtk")
