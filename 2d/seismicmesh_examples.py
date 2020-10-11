import SeismicMesh


def disk(h):
    bbox = (-1.0, 1.0, -1.0, 1.0)
    disk = SeismicMesh.geometry.Disk(0, 0, 1)

    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox,
        h0=h,
        domain=disk,
        edge_length=h,
        verbose=False
    )
    return points, cells


if __name__ == "__main__":
    import meshio

    points, cells = disk(0.1)
    meshio.Mesh(points, {"triangle": cells}).write("out.vtk")
