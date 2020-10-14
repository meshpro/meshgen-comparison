import numpy

import dmsh_examples
import meshpy_examples
import meshzoo_examples
import pygalmesh_examples
import pygmsh_examples
import seismicmesh_examples
from main import create_plots


def disk():
    functions = [
        pygalmesh_examples.disk,
        dmsh_examples.disk,
        pygmsh_examples.disk,
        meshpy_examples.disk,
        meshzoo_examples.disk,
        seismicmesh_examples.disk,
    ]
    # total runtime:
    # H = numpy.logspace(0.0, -1.5, num=15)  #  13.39s
    # H = numpy.logspace(0.0, -2.0, num=15)  # 227.95s
    # H = numpy.logspace(0.0, -2.1, num=15)  # 299.02s
    # H = numpy.logspace(-1.0, -2.1, num=15)  # 577.11s
    H = numpy.logspace(0.0, -2.1, num=15)  # 299.02s
    create_plots("disk", functions, H)


def l_shape():
    functions = [
        pygalmesh_examples.l_shape,
        dmsh_examples.l_shape,
        pygmsh_examples.l_shape,
        meshpy_examples.l_shape,
    ]
    H = numpy.logspace(0.0, -2.0, num=15)
    create_plots("l-shape", functions, H)


def rect_with_refinement():
    functions = [
        dmsh_examples.rect_with_refinement,
        pygmsh_examples.rect_with_refinement,
        meshpy_examples.rect_with_refinement,
        seismicmesh_examples.rect_with_refinement,
    ]
    # total runtime:
    H = numpy.logspace(0.0, -3.0, num=15)
    create_plots("rect-with-refinement", functions, H)


def ball():
    functions = [
        pygalmesh_examples.ball,
        pygmsh_examples.ball,
        meshpy_examples.ball,
        seismicmesh_examples.ball,
    ]
    # total runtime:
    H = numpy.logspace(0.0, -1.5, num=15)
    create_plots("ball", functions, H)


def l_shape_3d():
    functions = [pygalmesh_examples.l_shape_3d, pygmsh_examples.l_shape_3d]
    # total runtime:
    H = numpy.logspace(0.0, -1.5, num=15)
    create_plots("l-shape-3d", functions, H)


def box_with_refinement():
    functions = [
        pygalmesh_examples.box_with_refinement,
        pygmsh_examples.box_with_refinement,
        seismicmesh_examples.box_with_refinement,
    ]
    # total runtime:
    H = numpy.logspace(0.0, -2.0, num=15)
    create_plots("box-with-refinement", functions, H)


if __name__ == "__main__":
    disk()
    # l_shape()
    # rect_with_refinement()
    # ball()
    # box_with_refinement()
    # l_shape_3d()
    # l_shape()
