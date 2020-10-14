import dmsh
import meshpy
import meshzoo
import numpy
import pygalmesh
import pygmsh
import SeismicMesh

import dmsh_examples
import meshpy_examples
import meshzoo_examples
import pygalmesh_examples
import pygmsh_examples
import seismicmesh_examples
from main import create_plots


def disk():
    functions = {
        f"CGAL {pygalmesh.__cgal_version__}": pygalmesh_examples.disk,
        f"dmsh {dmsh.__version__}": dmsh_examples.disk,
        f"Gmsh {pygmsh.__gmsh_version__}": pygmsh_examples.disk,
        f"MeshPy {meshpy.version}": meshpy_examples.disk,
        f"meshzoo {meshzoo.__version__}": meshzoo_examples.disk,
        f"SeismicMesh {SeismicMesh.__version__}": seismicmesh_examples.disk,
    }
    # total runtime:
    # H = numpy.logspace(0.0, -1.5, num=15)  #  13.39s
    # H = numpy.logspace(0.0, -2.0, num=15)  # 227.95s
    # H = numpy.logspace(0.0, -2.1, num=15)  # 299.02s
    # H = numpy.logspace(-1.0, -2.1, num=15)  # 577.11s
    H = numpy.logspace(0.0, -2.1, num=15)  # 299.02s

    create_plots("disk", functions, H)


def l_shape():
    functions = {
        f"CGAL {pygalmesh.__cgal_version__}": pygalmesh_examples.l_shape,
        f"dmsh {dmsh.__version__}": dmsh_examples.l_shape,
        f"Gmsh {pygmsh.__gmsh_version__}": pygmsh_examples.l_shape,
        f"MeshPy {meshpy.version}": meshpy_examples.l_shape,
        # f"SeismicMesh {SeismicMesh.__version__}": seismicmesh_examples.disk,
    }
    # total runtime:
    H = numpy.logspace(0.0, -2.0, num=15)

    create_plots("l-shape", functions, H)


def rect_with_refinement():
    functions = {
        f"dmsh {dmsh.__version__}": dmsh_examples.rect_with_refinement,
        f"Gmsh {pygmsh.__gmsh_version__}": pygmsh_examples.rect_with_refinement,
        f"MeshPy {meshpy.version}": meshpy_examples.rect_with_refinement,
        f"SeismicMesh {SeismicMesh.__version__}": seismicmesh_examples.rect_with_refinement,
    }
    # total runtime:
    H = numpy.logspace(0.0, -3.0, num=15)

    create_plots("rect-with-refinement", functions, H)


def ball():
    functions = {
        f"CGAL {pygalmesh.__cgal_version__}": pygalmesh_examples.ball,
        f"Gmsh {pygmsh.__gmsh_version__}": pygmsh_examples.ball,
        f"MeshPy {meshpy.version}": meshpy_examples.ball,
        f"SeismicMesh {SeismicMesh.__version__}": seismicmesh_examples.ball,
    }
    # total runtime:
    H = numpy.logspace(0.0, -1.5, num=15)

    create_plots("ball", functions, H)


def box_with_refinement():
    functions = {
        f"CGAL {pygalmesh.__cgal_version__}": pygalmesh_examples.box_with_refinement,
        f"Gmsh {pygmsh.__gmsh_version__}": pygmsh_examples.box_with_refinement,
        # f"MeshPy {meshpy.version}": meshpy_examples.ball,
        f"SeismicMesh {SeismicMesh.__version__}": seismicmesh_examples.box_with_refinement,
    }
    # total runtime:
    H = numpy.logspace(0.0, -2.0, num=15)

    create_plots("box-with-refinement", functions, H)


if __name__ == "__main__":
    # disk()
    l_shape()
    # rect_with_refinement()
    # ball()
    # box_with_refinement()
