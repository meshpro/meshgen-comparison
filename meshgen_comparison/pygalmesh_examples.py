import numpy
import pygalmesh

desc = f"CGAL {pygalmesh.__cgal_version__} (via pygalmesh {pygalmesh.__version__})"
colors = ("#d62728", "#ff9896")


# CGAL may produce bad-quality meshes, see <https://github.com/CGAL/cgal/issues/5068>.
def disk(h):
    n = int(2 * numpy.pi / h)
    points = numpy.array(
        [
            [numpy.cos(alpha), numpy.sin(alpha)]
            for alpha in numpy.linspace(0.0, 2 * numpy.pi, n + 1, endpoint=False)
        ]
    )
    constraints = [[k, k + 1] for k in range(n)] + [[n, 0]]
    mesh = pygalmesh.generate_2d(
        points,
        constraints,
        # Relax max_edge size a bit; with 1.5, one gets node/cell numbers comparable to
        # the other mesh generators.
        max_edge_size=h * 1.5,
        num_lloyd_steps=0,
    )
    return mesh.points, mesh.get_cells_type("triangle")


def l_shape(h):
    points = numpy.array(
        [
            [-1.0, -1.0],
            [+1.0, -1.0],
            [+1.0, +0.0],
            [+0.0, +0.0],
            [+0.0, +1.0],
            [-1.0, +1.0],
        ]
    )
    n = len(points) - 1
    constraints = [[k, k + 1] for k in range(n)] + [[n, 0]]
    mesh = pygalmesh.generate_2d(
        points,
        constraints,
        # Relax max_edge size a bit; with 1.5, one gets node/cell numbers comparable to
        # the other mesh generators.
        max_edge_size=h * 1.5,
        num_lloyd_steps=0,
    )
    return mesh.points, mesh.get_cells_type("triangle")


def ball(h):
    s = pygalmesh.Ball([0, 0, 0], 1.0)
    # The circumradius of a regular tetrahedron with the given edge_size is sqrt(3 / 8)
    # * edge_size ~= 0.61 * edge_size. Relax it a bit and just use h.
    mesh = pygalmesh.generate_mesh(s, max_cell_circumradius=h, verbose=False)
    return mesh.points, mesh.get_cells_type("tetra")


def l_shape_3d(h):
    b0 = pygalmesh.Cuboid([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])
    b1 = pygalmesh.Cuboid([0.0, 0.0, 0.0], [1.1, 1.1, 1.1])
    u = pygalmesh.Difference(b0, b1)
    mesh = pygalmesh.generate_mesh(
        u,
        max_edge_size_at_feature_edges=h,
        max_cell_circumradius=h,
        verbose=False,
        # TODO sharpen the intersection edges and avoid
        # [1] 591479 segmentation fault (core dumped) python3 pygalmesh_examples.py
        # feature_edges=[
        #     [[1.0, 0.0, z] for z in numpy.linspace(0.0, 1.0, int(1.0 / h))]
        # ]
    )
    mesh.remove_lower_dimensional_cells()
    mesh.remove_orphaned_nodes()
    return mesh.points, mesh.get_cells_type("tetra")


def box_with_refinement(h):
    s = pygalmesh.Cuboid([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])

    def edge_length(x):
        return h + 0.1 * numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)

    mesh = pygalmesh.generate_mesh(
        s,
        max_edge_size_at_feature_edges=h,
        # The actual factor sqrt(3 / 8) leads to cells too small in comparison with
        # cells near the feature edges. Again, just take edge_size.
        max_cell_circumradius=edge_length,
        verbose=False,
    )

    return mesh.points, mesh.get_cells_type("tetra")


if __name__ == "__main__":
    import meshio

    points, cells = l_shape_3d(0.1)
    meshio.Mesh(points, {"tetra": cells}).write("out.vtk")
