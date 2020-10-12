import numpy
import pygalmesh


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
        # edge_size is the _max_ edge size. Relax it a bit; with 1.5, one gets node/cell
        # numbers comparable to the other mesh generators.
        edge_size=h * 1.5,
        num_lloyd_steps=0,
    )
    return mesh.points, mesh.get_cells_type("triangle")


def ball(h):
    s = pygalmesh.Ball([0, 0, 0], 1.0)
    # The circumradius of a regular tetrahedron with the given edge_size is
    # sqrt(3 / 8) * edge_size ~= 0.61 * edge_size. Relax it a bit and just use
    # edge_size.
    mesh = pygalmesh.generate_mesh(
        s,
        cell_size=h,
        verbose=False
    )
    return mesh.points, mesh.get_cells_type("tetra")


def box_with_refinement(h):
    s = pygalmesh.Cuboid([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])

    def edge_size(x):
        return h + 0.1 * numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)

    mesh = pygalmesh.generate_mesh(
        s,
        edge_size=edge_size,
        # The actual factor sqrt(3 / 8) leads to cells too small in comparison with
        # cells near the feature edges. Again, just take edge_size.
        cell_size=edge_size,
        verbose=False
    )

    return mesh.points, mesh.get_cells_type("tetra")


if __name__ == "__main__":
    import meshio

    points, cells = box_with_refinement(0.005)
    meshio.Mesh(points, {"tetra": cells}).write("out.vtk")
