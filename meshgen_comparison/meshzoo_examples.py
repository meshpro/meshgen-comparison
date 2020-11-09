import meshzoo

packages = [("meshzoo", meshzoo.__version__)]
colors = ("#8c564b", "#c49c94")  # cat20 brown


def disk(h):
    # tighten a bit
    h /= 1.1
    return meshzoo.disk(6, int(1 / h))


def sphere(h):
    # edge length of regular icosahedron with radius 1
    l = 1 / numpy.sin(0.4 * numpy.pi)
    n = int(l / h)
    return meshzoo.icosa_sphere(n)


if __name__ == "__main__":
    import numpy
    import meshio

    points, cells = sphere(0.3)
    meshio.Mesh(points, {"triangle": cells}).write("out.vtk")
