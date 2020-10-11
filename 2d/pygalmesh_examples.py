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


if __name__ == "__main__":
    disk(0.1)
