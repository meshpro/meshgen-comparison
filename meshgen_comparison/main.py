import datetime
import inspect
import json
import os
import signal
import time
import timeit
from contextlib import contextmanager
from pathlib import Path

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

from meshgen_comparison import (
    dmsh_examples,
    meshpy_examples,
    meshzoo_examples,
    pygalmesh_examples,
    pygmsh_examples,
    seismicmesh_examples,
)

# Some more colors:
# cat20_colors = [
#     ("#e377c2", "#f7b6d2"),
#     ("#7f7f7f", "#c7c7c7"),
#     ("#bcbd22", "#dbdb8d"),
#     ("#17becf", "#9edae5"),
# ]

modules = [
    dmsh_examples,
    meshpy_examples,
    meshzoo_examples,
    pygalmesh_examples,
    pygmsh_examples,
    seismicmesh_examples,
]
domains_h = [
    ("disk", numpy.logspace(0.0, -2.1, num=15)),  # 299.02s
    ("l_shape", numpy.logspace(0.0, -2.0, num=15)),
    ("rect_with_refinement", numpy.logspace(0.0, -3.0, num=15)),
    ("ball", numpy.logspace(-1.0, -3.0, num=15)),
    ("l_shape_3d", numpy.logspace(0.0, -1.5, num=15)),
    ("box_with_refinement", numpy.logspace(0.0, -2.0, num=15)),
]


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


def compute(fun, h):
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

    return {
        "time": toc - tic,
        # convert to python float types to JSON compatibility
        "quality_min": float(numpy.min(mesh.q_radius_ratio)),
        "quality_avg": float(numpy.average(mesh.q_radius_ratio)),
        "num_nodes": mesh.node_coords.shape[0],
        "num_poisson_steps": num_poisson_steps,
    }


def update_data_files():
    for module in modules:
        name = module.packages[0][0]
        print(name)

        # check if the existing data is up-to-date and skip if that's true
        json_filename = name.lower() + ".json"
        if Path(json_filename).is_file():
            with open(json_filename) as f:
                data = json.load(f)
            data_packages = [tuple(p) for p in data["packages"]]
            if set(module.packages) == set(data_packages):
                print(f"{json_filename} up to date.")
                print()
                continue

        # skip if computed before
        data = create_data(module)
        # store data in files
        with open(name.lower() + ".json", "w") as f:
            now = datetime.datetime.utcnow().replace(microsecond=0).isoformat()
            json.dump(
                {"name": name, "date": now, "packages": module.packages, "data": data},
                f,
                indent=2,
            )
        print()


def create_data(module, time_limit=120):
    # collect functions
    functions_h = []
    for domain, H in domains_h:
        try:
            fun = getattr(module, domain)
        except AttributeError:
            continue
        else:
            functions_h.append((fun, H))
    keys = [
        "h",
        "axpy_time",
        "time",
        "quality_min",
        "quality_avg",
        "num_nodes",
        "num_poisson_steps",
    ]
    data = {}
    with Progress() as progress:
        task0 = progress.add_task("domains", total=len(functions_h))
        task1 = progress.add_task("h")
        for fun, H in functions_h:
            data[fun.__name__] = {key: [] for key in keys}
            progress.update(task1, total=len(H))
            progress.reset(task1)
            for h in H:
                try:
                    with time_limiter(time_limit):
                        vals = compute(fun, h)
                except TimeoutException:
                    print("Timeout!", fun.__name__, h)
                    break
                for key, val in vals.items():
                    data[fun.__name__][key].append(val)
                data[fun.__name__]["h"].append(h)
                data[fun.__name__]["axpy_time"].append(_measure_axpy(vals["num_nodes"]))
                progress.update(task1, advance=1)
            progress.update(task0, advance=1)
    return data


def _measure_axpy(n):
    x = numpy.random.rand(n)
    y = numpy.random.rand(n)
    b = 3.14
    # get the minimum over 5 iterations
    return min(timeit.repeat(stmt=lambda: b * x + y, repeat=5, number=1))


def create_plots(domain):
    # plot the data
    plt.style.use(dufte.style)

    # num nodes vs. time
    for module in modules:
        name = module.packages[0][0]
        json_filename = name.lower() + ".json"
        with open(json_filename) as f:
            data = json.load(f)

        label = ", ".join(" ".join(package) for package in module.packages)
        if domain in data["data"]:
            d = data["data"][domain]
            plt.loglog(d["num_nodes"], d["time"], color=module.colors[0], label=label)

    dufte.legend()
    plt.xlabel("num points")
    plt.title("mesh creation times [s]")
    plt.savefig(f"{domain}-times.svg", transparent=True, bbox_inches="tight")
    # plt.show()
    plt.close()

    # num nodes vs. cell quality
    for module in modules:
        name = module.packages[0][0]
        json_filename = name.lower() + ".json"
        with open(json_filename) as f:
            data = json.load(f)
        label = ", ".join(" ".join(package) for package in module.packages)
        if domain in data["data"]:
            d = data["data"][domain]
            plt.semilogx(
                d["num_nodes"],
                d["quality_avg"],
                linestyle="-",
                color=module.colors[0],
                label=label,
            )
            plt.semilogx(
                d["num_nodes"], d["quality_min"], linestyle="--", color=module.colors[1]
            )
    dufte.legend()
    plt.xlabel("num points")
    plt.title("cell quality, avg  and min (dashed)")
    plt.ylim(0.0, 1.0)
    plt.savefig(f"{domain}-quality.svg", transparent=True, bbox_inches="tight")
    # plt.show()
    plt.close()

    # num nodes vs. CG iterations for Poisson
    for module in modules:
        name = module.packages[0][0]
        json_filename = name.lower() + ".json"
        with open(json_filename) as f:
            data = json.load(f)
        label = ", ".join(" ".join(package) for package in module.packages)
        if domain in data["data"]:
            d = data["data"][domain]
            plt.loglog(
                d["num_nodes"],
                d["num_poisson_steps"],
                color=module.colors[0],
                label=label,
            )
    dufte.legend()
    plt.xlabel("num points")
    poisson_tol = 1.0e-10
    plt.title(f"number of CG steps for the Poisson problem (tol={poisson_tol:.1e})")
    plt.savefig(f"{domain}-poisson.svg", transparent=True, bbox_inches="tight")
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
    update_data_files()
    for domain, _ in domains_h:
        create_plots(domain)
