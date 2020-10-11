import numpy
import SeismicMesh


def disk(h):
    bbox = (-1.0, 1.0, -1.0, 1.0)
    domain = SeismicMesh.geometry.Disk(0, 0, 1)

    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox, h0=h, domain=domain, edge_length=h
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


if __name__ == "__main__":
    import meshio

    # points, cells = rect_with_refinement(0.1)
    points, cells = ball(0.1)
    meshio.Mesh(points, {"tetra": cells}).write("out.vtk")