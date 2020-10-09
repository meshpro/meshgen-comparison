import pygmsh
import math


def disk(h):
    with pygmsh.geo.Geometry() as geom:
        geom.add_circle(
            [0.0, 0.0, 0.0],
            1.0,
            mesh_size=h,
            num_sections=4,
            compound=True,
        )
        mesh = geom.generate_mesh()
    mesh.prune_z_0()
    return mesh.points, mesh.get_cells_type("triangle")


def rect_with_refinement(h):
    with pygmsh.geo.Geometry() as geom:
        geom.add_polygon([
            [-1.0, -1.0],
            [+1.0, -1.0],
            [+1.0, +1.0],
            [-1.0, +1.0],
        ])
        geom.set_mesh_size_callback(
            lambda dim, tag, x, y, z: h + 0.1 * math.sqrt(x ** 2 + y ** 2)
        )
        mesh = geom.generate_mesh()

    mesh.prune_z_0()
    return mesh.points, mesh.get_cells_type("triangle")


if __name__ == "__main__":
    import meshio

    points, cells = rect_with_refinement(0.01)
    meshio.Mesh(points, {"triangle": cells}).write("out.vtk")
