import math

import meshio
import meshpy.triangle
import numpy

packages = [("MeshPy", meshpy.version)]
colors = ("#1f77b4", "#aec7e8")  # cat20 blue


def _round_trip_connect(end):
    return [[i, i + 1] for i in range(end)] + [[end, 0]]


def disk(h):
    num_boundary_points = int(2 * numpy.pi / h)
    n_phi = num_boundary_points
    radius = 1.0

    # Choose the maximum area of a triangle equal to the area of an equilateral triangle
    # on the boundary.
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


def l_shape(h):
    points = [
        [-1.0, -1.0],
        [+1.0, -1.0],
        [+1.0, +0.0],
        [+0.0, +0.0],
        [+0.0, +1.0],
        [-1.0, +1.0],
    ]

    info = meshpy.triangle.MeshInfo()
    info.set_points(points)
    info.set_facets(_round_trip_connect(len(points) - 1))

    # Choose the maximum area of a triangle equal to the area of an equilateral triangle
    # with the edge length h.
    max_area = math.sqrt(3) / 4 * h ** 2
    # relax a bit
    max_area *= 1.5

    mesh = meshpy.triangle.build(info, refinement_func=lambda v, a: a > max_area)

    return numpy.array(mesh.points), numpy.array(mesh.elements)


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
        lim = h + 0.1 * math.sqrt(bc[0] ** 2 + bc[1] ** 2)
        # relax the limit a bit; the refinement is too strict otherwise
        return numpy.any(edge_lengths > lim * 1.5)

    info = meshpy.triangle.MeshInfo()
    info.set_points(points)
    info.set_facets(_round_trip_connect(len(points) - 1))

    mesh = meshpy.triangle.build(info, refinement_func=needs_refinement)

    return numpy.array(mesh.points), numpy.array(mesh.elements)


# segfaults, doesn't react to term signals
# https://github.com/inducer/meshpy/issues/58
# def ball(h):
#     from meshpy.geometry import (
#         EXT_OPEN,
#         GeometryBuilder,
#         generate_surface_of_revolution,
#     )
#     from meshpy.tet import MeshInfo, build
#
#     r = 3
#
#     polar_subdivision = int(math.pi / h)
#     dphi = math.pi / polar_subdivision
#
#     def truncate(val):
#         return 0 if abs(val) < 1e-10 else val
#
#     rz = [
#         [truncate(r * math.sin(i * dphi)), r * math.cos(i * dphi)]
#         for i in range(polar_subdivision + 1)
#     ]
#
#     geob = GeometryBuilder()
#     radial_subdivision = int(2 * math.pi / h)
#     geob.add_geometry(
#         *generate_surface_of_revolution(
#             rz, closure=EXT_OPEN, radial_subdiv=radial_subdivision
#         )
#     )
#
#     mesh_info = MeshInfo()
#     geob.set(mesh_info)
#
#     mesh = build(mesh_info)
#
#     # remove orphaned nodes
#     mesh = meshio.Mesh(numpy.array(mesh.points), {"tetra": numpy.array(mesh.elements)})
#     mesh.remove_orphaned_nodes()
#
#     return mesh.points, mesh.get_cells_type("tetra")


if __name__ == "__main__":
    points, cells = l_shape(0.1)
    # points, cells = ball(0.1)
    meshio.Mesh(points, {"triangle": cells}).write("out.vtk")
