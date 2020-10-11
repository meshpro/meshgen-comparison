import numpy
import SeismicMesh


def disk(h):
    bbox = (-1.0, 1.0, -1.0, 1.0)
    disk = SeismicMesh.geometry.Disk(0, 0, 1)

    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox, h0=h, domain=disk, edge_length=h, verbose=False
    )
    return points, cells


def rect_with_refinement(h):
    bbox = (-1.0, 1.0, -1.0, 1.0)
    rect = SeismicMesh.geometry.Rectangle([-1.0, +1.0, -1.0, +1.0])

    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox,
        h0=h,
        domain=rect,
        edge_length=lambda x: h + 0.1 * numpy.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2),
        verbose=False,
    )
    return points, cells


if __name__ == "__main__":
    import meshio

    points, cells = rect_with_refinement(0.1)
    meshio.Mesh(points, {"triangle": cells}).write("out.vtk")
