import math

import meshpy.triangle
import numpy


def _round_trip_connect(end):
    return [[i, i + 1] for i in range(end)] + [[end, 0]]


def disk(h):
    num_boundary_points = int(2 * numpy.pi / h)
    n_phi = num_boundary_points
    radius = 1.0

    # Choose the maximum area of a triangle equal to the area of
    # an equilateral triangle on the boundary.
    a_boundary = 2 * numpy.pi * radius / n_phi
    max_area = a_boundary ** 2 * numpy.sqrt(3.0) / 4.0
    max_area = float(max_area)  # meshpy can't deal with numpy.float64
    # relax a bit
    max_area *= 1.5

    # generate points on the boundary
    boundary_points = [
        [radius * numpy.cos(phi), radius * numpy.sin(phi)]
        for phi in numpy.linspace(0.0, 2 * numpy.pi, n_phi, endpoint=False)
    ]

    # create the mesh
    info = meshpy.triangle.MeshInfo()
    info.set_points(boundary_points)

    info.set_facets(_round_trip_connect(len(boundary_points) - 1))

    def _needs_refinement(vertices, area):
        return area > max_area

    meshpy_mesh = meshpy.triangle.build(info, refinement_func=_needs_refinement)

    return numpy.array(meshpy_mesh.points), numpy.array(meshpy_mesh.elements)


def rect_with_refinement(h):
    points = [(1, 1), (-1, 1), (-1, -1), (1, -1)]

    def needs_refinement(verts, area):
        edges = numpy.array(
            [
                [verts[0].x - verts[1].x, verts[0].y - verts[1].y],
                [verts[1].x - verts[2].x, verts[1].y - verts[2].y],
                [verts[2].x - verts[0].x, verts[2].y - verts[0].y],
            ]
        )
        edge_lengths = numpy.sqrt(numpy.einsum("ij,ij->i", edges, edges))
        bc = [
            (verts[0].x + verts[1].x + verts[2].x) / 3,
            (verts[0].y + verts[1].y + verts[2].y) / 3,
        ]
        return numpy.any(edge_lengths > h + 0.1 * math.sqrt(bc[0] ** 2 + bc[1] ** 2))

    info = meshpy.triangle.MeshInfo()
    info.set_points(points)
    info.set_facets(_round_trip_connect(len(points) - 1))

    mesh = meshpy.triangle.build(info, refinement_func=needs_refinement)

    points = numpy.array(mesh.points)
    cells = numpy.array(mesh.elements)

    # import meshio
    # meshio.Mesh(points, {"triangle": cells}).write("out.vtk")

    return points, cells


if __name__ == "__main__":
    rect_with_refinement(0.1)
