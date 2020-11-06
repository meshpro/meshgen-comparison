import math

import numpy
import pygmsh

packages = [
    ("Gmsh", pygmsh.__gmsh_version__),
    ("pygmsh", pygmsh.__version__),
]
colors = ("#d62728", "#ff9896")  # cat20 red


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


def l_shape(h):
    with pygmsh.geo.Geometry() as geom:
        geom.add_polygon(
            [
                [-1.0, -1.0],
                [+1.0, -1.0],
                [+1.0, +0.0],
                [+0.0, +0.0],
                [+0.0, +1.0],
                [-1.0, +1.0],
            ],
            mesh_size=h,
        )
        mesh = geom.generate_mesh()

    mesh.prune_z_0()
    return mesh.points, mesh.get_cells_type("triangle")


def rect_with_refinement(h):
    with pygmsh.geo.Geometry() as geom:
        geom.add_polygon(
            [
                [-1.0, -1.0],
                [+1.0, -1.0],
                [+1.0, +1.0],
                [-1.0, +1.0],
            ]
        )
        geom.set_mesh_size_callback(
            lambda dim, tag, x, y, z: h + 0.1 * math.sqrt(x ** 2 + y ** 2)
        )
        mesh = geom.generate_mesh()

    mesh.prune_z_0()
    return mesh.points, mesh.get_cells_type("triangle")


def quarter_annulus(h):
    with pygmsh.occ.Geometry() as geom:
        r0 = 0.25
        r1 = 1.0
        disk0 = geom.add_disk([0.0, 0.0], r0)
        disk1 = geom.add_disk([0.0, 0.0], r1)
        diff0 = geom.boolean_difference(disk1, disk0)

        rect = geom.add_rectangle([0.0, 0.0, 0.0], 1.0, 1.0)
        geom.boolean_intersection([diff0, rect])

        geom.set_mesh_size_callback(
            lambda dim, tag, x, y, z: h + 0.1 * (numpy.sqrt(x ** 2 + y ** 2) - r0),
        )
        mesh = geom.generate_mesh()

    mesh.prune_z_0()
    return mesh.points, mesh.get_cells_type("triangle")


def ball(h):
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_max = h
        geom.add_ball([0, 0, 0], 1.0)
        mesh = geom.generate_mesh()
    return mesh.points, mesh.get_cells_type("tetra")


def cylinder(h):
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_max = h
        geom.add_cylinder([0, 0, 0], [0.0, 0.0, 1.0], 0.5)
        mesh = geom.generate_mesh()
    return mesh.points, mesh.get_cells_type("tetra")


def l_shape_3d(h):
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_max = h
        b0 = geom.add_box([-1.0, -1.0, -1.0], [2.0, 2.0, 2.0])
        b1 = geom.add_box([0.0, 0.0, 0.0], [2.0, 2.0, 2.0])
        geom.boolean_difference(b0, b1)
        mesh = geom.generate_mesh()
    # mesh.remove_orphaned_nodes()
    return mesh.points, mesh.get_cells_type("tetra")


def box_with_refinement(h):
    with pygmsh.geo.Geometry() as geom:
        geom.add_box(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)
        geom.set_mesh_size_callback(
            lambda dim, tag, x, y, z: h + 0.1 * math.sqrt(x ** 2 + y ** 2 + z ** 2)
        )
        mesh = geom.generate_mesh()

    return mesh.points, mesh.get_cells_type("tetra")


if __name__ == "__main__":
    import meshio

    # points, cells = l_shape(0.1)
    # points, cells = rect_with_refinement(0.01)
    # points, cells = l_shape_3d(0.1)
    # points, cells = box_with_refinement(0.01)
    # points, cells = cylinder(0.05)
    points, cells = quarter_annulus(0.01)
    meshio.Mesh(points, {"triangle": cells}).write("out.vtk")
