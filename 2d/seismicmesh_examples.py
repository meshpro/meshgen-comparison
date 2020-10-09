import SeismicMesh
import numpy


def disk(h):
    bbox = (-1.0, 1.0, -1.0, 1.0)
    circle = SeismicMesh.geometry.Circle(0, 0, 1)

    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox,
        h0=h,
        domain=circle,
        cell_size=lambda x: numpy.full(len(x), h),
        # verbose=False
    )
    return points, cells


if __name__ == "__main__":
    import meshio

    points, cells = disk(0.1)
    meshio.Mesh(points, {"triangle": cells}).write("out.vtk")
