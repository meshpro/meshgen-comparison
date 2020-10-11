import os
import time

import dmsh
import dufte
import matplotlib.pyplot as plt
import meshio
import meshplex
import meshpy
import meshzoo
import numpy
import pygalmesh
import pygmsh
import SeismicMesh
from dolfin import (
    Constant,
    DirichletBC,
    Function,
    FunctionSpace,
    KrylovSolver,
    Mesh,
    TestFunction,
    TrialFunction,
    XDMFFile,
    assemble,
    dx,
    grad,
    inner,
)
from rich.progress import Progress

import dmsh_examples
import meshpy_examples
import meshzoo_examples
import pygalmesh_examples
import pygmsh_examples
import seismicmesh_examples

cat20_colors = [
    ("#1f77b4", "#aec7e8"),
    ("#ff7f0e", "#ffbb78"),
    ("#2ca02c", "#98df8a"),
    ("#d62728", "#ff9896"),
    ("#9467bd", "#c5b0d5"),
    ("#8c564b", "#c49c94"),
    ("#e377c2", "#f7b6d2"),
    ("#7f7f7f", "#c7c7c7"),
    ("#bcbd22", "#dbdb8d"),
    ("#17becf", "#9edae5"),
]


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


def rect_with_refinement():
    functions = {
        f"dmsh {dmsh.__version__}": dmsh_examples.rect_with_refinement,
        f"Gmsh {pygmsh.__gmsh_version__}": pygmsh_examples.rect_with_refinement,
        f"MeshPy {meshpy.version}": meshpy_examples.rect_with_refinement,
        f"SeismicMesh {SeismicMesh.__version__}": seismicmesh_examples.rect_with_refinement,
    }
    # total runtime:
    H = numpy.logspace(0.0, -2.1, num=15)

    create_plots("rect-with-refinement", functions, H)


def ball():
    functions = {
        f"CGAL {pygalmesh.__cgal_version__}": pygalmesh_examples.ball,
        f"Gmsh {pygmsh.__gmsh_version__}": pygmsh_examples.ball,
        # f"MeshPy {meshpy.version}": meshpy_examples.rect_with_refinement,
        f"SeismicMesh {SeismicMesh.__version__}": seismicmesh_examples.ball,
    }
    # total runtime:
    H = numpy.logspace(0.0, -2.1, num=15)

    create_plots("ball", functions, H)


def create_plots(prefix, functions, H):
    times = []
    quality_min = []
    quality_avg = []
    num_poisson_steps = []
    num_points = []

    colors = cat20_colors[: len(functions)]

    poisson_tol = 1.0e-10
    with Progress() as progress:
        task1 = progress.add_task("Overall", total=len(H))
        task2 = progress.add_task("Functions", total=len(functions))
        for h in H:
            times.append([])
            quality_min.append([])
            quality_avg.append([])
            num_poisson_steps.append([])
            num_points.append([])
            progress.update(task2, completed=0)
            for fun in functions.values():
                tic = time.time()
                points, cells = fun(h)
                toc = time.time()

                if cells.shape[1] == 3:
                    mesh = meshplex.MeshTri(points, cells)
                else:
                    assert cells.shape[1] == 4
                    mesh = meshplex.MeshTetra(points, cells)

                times[-1].append(toc - tic)
                quality_min[-1].append(numpy.min(mesh.q_radius_ratio))
                quality_avg[-1].append(numpy.average(mesh.q_radius_ratio))
                num_points[-1].append(mesh.node_coords.shape[0])

                if numpy.min(mesh.q_radius_ratio) < 1.0e-5:
                    num_poisson_steps[-1].append(numpy.nan)
                else:
                    num_steps = get_poisson_steps(points, cells, poisson_tol)
                    num_poisson_steps[-1].append(num_steps)

                progress.update(task2, advance=1)
            progress.update(task1, advance=1)

    times = numpy.array(times)
    quality_min = numpy.array(quality_min)
    quality_avg = numpy.array(quality_avg)
    num_poisson_steps = numpy.array(num_poisson_steps)
    num_points = numpy.array(num_points)

    # plot the data
    plt.style.use(dufte.style)
    for name, num_pts, t, cols in zip(functions.keys(), num_points.T, times.T, colors):
        plt.loglog(num_pts, t, color=cols[0], label=name)
    dufte.legend()
    plt.xlabel("num points")
    plt.title("mesh creation times [s]")
    plt.savefig(f"{prefix}-times.svg", transparent=True, bbox_inches="tight")
    # plt.show()
    plt.close()

    for name, num_pts, qa, qm, cols in zip(
        functions.keys(), num_points.T, quality_avg.T, quality_min.T, colors
    ):
        plt.semilogx(num_pts, qa, color=cols[0], linestyle="-", label=f"{name}")
        plt.semilogx(num_pts, qm, color=cols[1], linestyle="--")
        plt.ylim(0.0, 1.0)
    dufte.legend()
    plt.xlabel("num points")
    plt.title("cell quality, avg  and min (dashed)")
    plt.savefig(f"{prefix}-quality.svg", transparent=True, bbox_inches="tight")
    # plt.show()
    plt.close()

    for name, num_pts, np, cols in zip(
        functions.keys(), num_points.T, num_poisson_steps.T, colors
    ):
        plt.semilogx(num_pts, np, color=cols[0], label=f"{name}")
    dufte.legend()
    plt.xlabel("num points")
    plt.title(f"number of CG steps for the Poisson problem (tol={poisson_tol:.1e})")
    plt.savefig(f"{prefix}-poisson.svg", transparent=True, bbox_inches="tight")
    # plt.show()
    plt.close()


def get_poisson_steps(pts, cells, tol):
    # Still can't initialize a mesh from points, cells
    filename = "mesh.xdmf"
    if cells.shape[1] == 3:
        meshio.write_points_cells(filename, pts, {"triangle": cells})
    else:
        assert cells.shape[1] == 4
        meshio.write_points_cells(filename, pts, {"tetra": cells})

    mesh = Mesh()
    with XDMFFile(filename) as f:
        f.read(mesh)
    os.remove(filename)
    os.remove("mesh.h5")

    # build Laplace matrix with Dirichlet boundary using dolfin
    V = FunctionSpace(mesh, "Lagrange", 1)

    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(grad(u), grad(v)) * dx
    u0 = Constant(0.0)
    bc = DirichletBC(V, u0, "on_boundary")
    f = Constant(1.0)
    L = f * v * dx

    A = assemble(a)
    b = assemble(L)
    bc.apply(A, b)

    # solve(A, x, b, "cg")
    solver = KrylovSolver("cg", "none")
    solver.parameters["absolute_tolerance"] = 0.0
    solver.parameters["relative_tolerance"] = tol

    x = Function(V)
    x_vec = x.vector()
    num_steps = solver.solve(A, x_vec, b)

    # # convert to scipy matrix
    # A = as_backend_type(A).mat()
    # ai, aj, av = A.getValuesCSR()
    # A = scipy.sparse.csr_matrix(
    #     (av, aj, ai), shape=(A.getLocalSize()[0], A.getSize()[1])
    # )

    # # ev = eigvals(A.todense())
    # ev_max = scipy.sparse.linalg.eigs(A, k=1, which="LM")[0][0]
    # assert numpy.abs(ev_max.imag) < 1.0e-15
    # ev_max = ev_max.real
    # ev_min = scipy.sparse.linalg.eigs(A, k=1, which="SM")[0][0]
    # assert numpy.abs(ev_min.imag) < 1.0e-15
    # ev_min = ev_min.real
    # cond = ev_max / ev_min

    # solve poisson system, count num steps
    # b = numpy.ones(A.shape[0])
    # out = pykry.gmres(A, b)
    # num_steps = len(out.resnorms)
    return num_steps


if __name__ == "__main__":
    # disk()
    # rect_with_refinement()
    ball()
