import inspect
import json
import os
import signal
import time
from contextlib import contextmanager

import dufte
import matplotlib.pyplot as plt
import meshio
import meshplex
import numpy
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

# Some more colors:
# cat20_colors = [
#     ("#e377c2", "#f7b6d2"),
#     ("#7f7f7f", "#c7c7c7"),
#     ("#bcbd22", "#dbdb8d"),
#     ("#17becf", "#9edae5"),
# ]


# https://stackoverflow.com/a/601168/353337
class TimeoutException(Exception):
    pass


@contextmanager
def time_limiter(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")

    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


def compute(name, fun, h):
    tic = time.time()
    points, cells = fun(h)
    toc = time.time()

    if cells.shape[1] == 3:
        mesh = meshplex.MeshTri(points, cells)
    else:
        assert cells.shape[1] == 4
        mesh = meshplex.MeshTetra(points, cells)

    if numpy.min(mesh.q_radius_ratio) < 1.0e-5:
        num_poisson_steps = numpy.nan
    else:
        poisson_tol = 1.0e-10
        num_poisson_steps = get_poisson_steps(points, cells, poisson_tol)

    return (
        toc - tic,
        numpy.min(mesh.q_radius_ratio),
        numpy.average(mesh.q_radius_ratio),
        mesh.node_coords.shape[0],
        num_poisson_steps,
    )


def create_data(domain, functions, H, time_limit=60):
    names = [inspect.getmodule(fun).packages[0][0].tolower() for fun in functions]
    keys = ["h", "time", "quality_min", "quality_avg", "num_nodes", "num_poisson_steps"]
    data = {name: {key: [] for key in keys} for name in names}

    with Progress() as progress:
        task1 = progress.add_task("Functions", total=len(functions))
        task2 = progress.add_task("h", total=len(H))
        for name, fun in zip(names, functions):
            progress.update(task2, completed=0)
            for h in H:
                try:
                    with time_limiter(time_limit):
                        vals = compute(name, fun, h)
                except TimeoutException:
                    print("Timeout!", name, h)
                    break
                assert len(keys[1:]) == len(vals)
                for key, val in zip(keys[1:], vals):
                    data[name][key].append(val)
                data[name]["h"].append(h)
                progress.update(task2, advance=1)
            progress.update(task1, advance=1)

    # store data in files
    with open(domain + ".json", "w") as f:
        json.dump(data, f, indent=2)
    exit(1)


def create_plots():
    colors = [inspect.getmodule(fun).colors for fun in functions]

    # plot the data
    plt.style.use(dufte.style)
    for name, num_pts, t, cols in zip(names, num_points.T, times.T, colors):
        plt.loglog(num_pts, t, color=cols[0], label=name)
    dufte.legend()
    plt.xlabel("num points")
    plt.title("mesh creation times [s]")
    plt.savefig(f"{prefix}-times.svg", transparent=True, bbox_inches="tight")
    # plt.show()
    plt.close()

    for name, num_pts, qa, qm, cols in zip(
        names, num_points.T, quality_avg.T, quality_min.T, colors
    ):
        plt.semilogx(num_pts, qa, color=cols[0], linestyle="-", label=f"{name}")
        plt.semilogx(num_pts, qm, color=cols[1], linestyle="--", label="")
        plt.ylim(0.0, 1.0)
    dufte.legend()
    plt.xlabel("num points")
    plt.title("cell quality, avg  and min (dashed)")
    plt.savefig(f"{prefix}-quality.svg", transparent=True, bbox_inches="tight")
    # plt.show()
    plt.close()

    for name, num_pts, np, cols in zip(
        names, num_points.T, num_poisson_steps.T, colors
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
