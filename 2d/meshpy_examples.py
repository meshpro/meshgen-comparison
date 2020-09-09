import meshpy.triangle
import numpy


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
    Phi = numpy.linspace(0.0, 2 * numpy.pi, n_phi, endpoint=False)
    boundary_points = [
        [radius * numpy.cos(phi), radius * numpy.sin(phi)] for phi in Phi
    ]

    # create the mesh
    info = meshpy.triangle.MeshInfo()
    info.set_points(boundary_points)

    def _round_trip_connect(end):
        return [[i, i + 1] for i in range(end)] + [[end, 0]]

    info.set_facets(_round_trip_connect(len(boundary_points) - 1))

    def _needs_refinement(vertices, area):
        return area > max_area

    meshpy_mesh = meshpy.triangle.build(info, refinement_func=_needs_refinement)

    return numpy.array(meshpy_mesh.points), numpy.array(meshpy_mesh.elements)
