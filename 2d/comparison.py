import time

from dolfin import (
    Mesh,
    XDMFFile,
    TrialFunction,
    TestFunction,
    Function,
    FunctionSpace,
    inner,
    grad,
    Constant,
    dx,
    DirichletBC,
    assemble,
    KrylovSolver,
)
import dufte
import meshio
import meshplex
import matplotlib.pyplot as plt
import numpy
from rich.progress import Progress

import dmsh_examples
import meshpy_examples
import meshzoo_examples
import pygalmesh_examples
import pygmsh_examples


def disk():
    # total runtime:
    # H = numpy.logspace(0, -1.5, num=15)  #  13.39s
    # H = numpy.logspace(0, -2.0, num=15)  # 227.95s
    # H = numpy.logspace(0, -2.1, num=15)  # 299.02s
    # H = numpy.logspace(-1.0, -2.1, num=15)  # 577.11s
    H = numpy.logspace(0, -2.1, num=15)  # 299.02s
    times = []
    quality_min = []
    quality_avg = []
    num_poisson_steps = []
    num_points = []
    modules = {
        "cgal": pygalmesh_examples,
        "dmsh": dmsh_examples,
        "gmsh": pygmsh_examples,
        "meshpy": meshpy_examples,
        "meshzoo": meshzoo_examples,
    }
    # cat20 colors
    colors = [
        ("#1f77b4", "#aec7e8"),
        ("#ff7f0e", "#ffbb78"),
        ("#2ca02c", "#98df8a"),
        ("#d62728", "#ff9896"),
        ("#9467bd", "#c5b0d5"),
        # ("#8c564b", "#c49c94"),
        # ("#e377c2", "#f7b6d2"),
        # ("#7f7f7f", "#c7c7c7"),
        # ("#bcbd22", "#dbdb8d"),
        # ("#17becf", "#9edae5"),
    ]
    assert len(colors) == len(modules)

    poisson_tol = 1.0e-10
    with Progress() as progress:
        task1 = progress.add_task("Overall", total=len(H))
        task2 = progress.add_task("Modules", total=len(modules))
        for h in H:
            times.append([])
            quality_min.append([])
            quality_avg.append([])
            num_poisson_steps.append([])
            num_points.append([])
            progress.update(task2, completed=0)
            for module in modules.values():
                tic = time.time()
                points, cells = module.disk(h)
                toc = time.time()

                mesh = meshplex.MeshTri(points, cells)

                times[-1].append(toc - tic)
                quality_min[-1].append(numpy.min(mesh.cell_quality))
                quality_avg[-1].append(numpy.average(mesh.cell_quality))
                num_points[-1].append(mesh.node_coords.shape[0])

                if numpy.min(mesh.cell_quality) < 1.0e-5:
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

    print("names:")
    print(list(modules.keys()))
    print()
    print("num_points:")
    print(num_points)
    print()
    print("times:")
    print(times)
    print()
    print("quality_avg:")
    print(quality_avg)
    print()
    print("quality_min:")
    print(quality_min)
    print()
    print("num_poisson_steps:")
    print(num_poisson_steps)

    # plot the data
    plt.style.use(dufte.style)
    for name, num_pts, t, cols in zip(modules.keys(), num_points.T, times.T, colors):
        plt.loglog(num_pts, t, color=cols[0], label=name)
    dufte.legend()
    plt.xlabel("num points")
    plt.title("mesh creation times [s]")
    plt.savefig("disk-times.svg", transparent=True, bbox_inches="tight")
    # plt.show()
    plt.close()

    for name, num_pts, qa, qm, cols in zip(
        modules.keys(), num_points.T, quality_avg.T, quality_min.T, colors
    ):
        plt.semilogx(num_pts, qa, color=cols[0], linestyle="-", label=f"{name} (avg)")
        plt.semilogx(num_pts, qm, color=cols[1], linestyle="--", label=f"{name} (min)")
        plt.ylim(0.0, 1.0)
    dufte.legend()
    plt.xlabel("num points")
    plt.title("cell quality")
    plt.savefig("disk-quality.svg", transparent=True, bbox_inches="tight")
    # plt.show()
    plt.close()

    for name, num_pts, np, cols in zip(
        modules.keys(), num_points.T, num_poisson_steps.T, colors
    ):
        plt.semilogx(num_pts, np, color=cols[0], label=f"{name}")
    dufte.legend()
    plt.xlabel("num points")
    plt.title(f"number of CG steps for the Poisson problem (tol={poisson_tol:.1e})")
    plt.savefig("disk-poisson.svg", transparent=True, bbox_inches="tight")
    # plt.show()
    plt.close()


def get_poisson_steps(pts, cells, tol):
    # Still can't initialize a mesh from points, cells
    filename = "mesh.xdmf"
    meshio.write_points_cells(filename, pts, {"triangle": cells})
    mesh = Mesh()
    with XDMFFile(filename) as f:
        f.read(mesh)

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
    solver.parameters['absolute_tolerance'] = 0.0
    solver.parameters['relative_tolerance'] = tol

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
    disk()
