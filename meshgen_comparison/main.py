import datetime
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
import quadpy
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
# colors = ("#1f77b4", "#aec7e8")  # cat20 blue
# colors = ("#ff7f0e", "#ffbb78")  # cat20 orange
# colors = ("#2ca02c", "#98df8a")  # cat20 green
# colors = ("#d62728", "#ff9896")  # cat20 red
# colors = ("#9467bd", "#c5b0d5")  # cat20 purple
# colors = ("#8c564b", "#c49c94")  # cat20 brown
#     ("#e377c2", "#f7b6d2"),  # pink
#     ("#7f7f7f", "#c7c7c7"),  # gray
#     ("#bcbd22", "#dbdb8d"),  # yellow
#     ("#17becf", "#9edae5"),  # light blue  # gmsh?
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
    ("disk", numpy.logspace(-1.0, -2.1, num=15)),
    ("l_shape", numpy.logspace(-1.0, -2.1, num=15)),
    ("rect_with_refinement", numpy.logspace(-1.0, -3.0, num=15)),
    ("quarter_annulus", numpy.logspace(-1.0, -3.0, num=15)),
    ("sphere", numpy.logspace(-1.0, -2.5, num=15)),
    ("ball", numpy.logspace(-1.0, -3.0, num=15)),
    ("cylinder", numpy.logspace(-1.0, -1.8, num=15)),
    ("l_shape_3d", numpy.logspace(-1.0, -1.4, num=15)),
    ("box_with_refinement", numpy.logspace(-1.0, -2.0, num=15)),
]


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


def energy(mesh):
    """This is the mesh energy given by

       E = int_Omega |u_l - u| * rho

    with u = ||x||^2 and u_l being its piecewise linearization. This energy (or an
    approximation thereof) is minimized in mesh optimization.

    Note that u_l >= u, so the absolute value can be removed.

    Only works for rho=1 right now.
    """
    if isinstance(mesh, meshplex.MeshTri):
        scheme = quadpy.t2.get_good_scheme(2)
    else:
        assert isinstance(mesh, meshplex.MeshTetra)
        scheme = quadpy.t3.get_good_scheme(2)

    triangles = numpy.moveaxis(mesh.points[mesh.cells["points"]], 0, 1)

    def u(x):
        return numpy.einsum("i...,i...->...", x, x)

    # int_Omega u
    vals = scheme.integrate(u, triangles)
    int_u = numpy.sum(vals)
    # if uniform_density:
    #     val = numpy.sum(val)
    # else:
    #     rho = 1.0 / mesh.cell_volumes
    #     val = numpy.dot(val, rho)

    # vertex scheme int_Omega u_l
    vals = numpy.einsum("ij,ij->i", mesh.points, mesh.points)
    int_ul = numpy.sum(
        numpy.sum(vals[mesh.cells["points"]], axis=1) / 3 * mesh.cell_volumes
    )
    return int_ul - int_u


def compute(fun, h):
    tic = time.time()
    points, cells = fun(h)
    toc = time.time()

    if cells.shape[1] == 3:
        mesh = meshplex.MeshTri(points, cells)
    else:
        assert cells.shape[1] == 4
        mesh = meshplex.MeshTetra(points, cells)

    if fun.__name__ == "sphere" or numpy.min(mesh.q_radius_ratio) < 1.0e-5:
        num_poisson_steps = numpy.nan
    else:
        poisson_tol = 1.0e-10
        # num_poisson_steps = get_poisson_steps_dolfin(points, cells, poisson_tol)
        num_poisson_steps = get_poisson_steps_scikitfem(points, cells, poisson_tol)

    return {
        "time": toc - tic,
        # convert to python float types to JSON compatibility
        "quality_min": float(numpy.min(mesh.q_radius_ratio)),
        "quality_avg": float(numpy.average(mesh.q_radius_ratio)),
        "energy": energy(mesh),
        "num_nodes": mesh.points.shape[0],
        "num_poisson_steps": num_poisson_steps,
    }


def update_data_file(module):
    module_name, version = module.packages[0]
    print(f"{module_name} {version}")

    # collect all possible functions
    functions_h = []
    for domain, H in domains_h:
        try:
            fun = getattr(module, domain)
        except AttributeError:
            continue
        else:
            functions_h.append((fun, H))

    # check if the existing data is up-to-date and skip if that's true
    json_filename = module_name.lower() + ".json"
    if Path(json_filename).is_file():
        with open(json_filename) as f:
            content = json.load(f)
        data_packages = [tuple(p) for p in content["packages"]]
        if set(module.packages) == set(data_packages):
            # versions up-to-date; check if all domains are present
            functions_h = [
                fun_h
                for fun_h in functions_h
                if fun_h[0].__name__ not in content["data"]
            ]
        else:
            # version outdated; rerun all
            content = {"data": {}}
    else:
        content = {"data": {}}

    if len(functions_h) == 0:
        print(f"{json_filename} up to date.")
        print()
        return

    string = ", ".join(f.__name__ for f, _ in functions_h)
    print(f"{len(functions_h)} domains remaining ({string}):")

    keys = [
        "h",
        "axpy_time",
        "time",
        "quality_min",
        "quality_avg",
        "energy",
        "num_nodes",
        "num_poisson_steps",
    ]
    time_limit = 120
    with Progress() as progress:
        task0 = progress.add_task("domains", total=len(functions_h))
        task1 = progress.add_task("h")
        for fun, H in functions_h:
            data = {key: [] for key in keys}
            progress.update(task1, total=len(H))
            progress.reset(task1)
            for h in H:
                try:
                    with time_limiter(time_limit):
                        vals = compute(fun, h)
                except TimeoutException:
                    print("Timeout!", fun.__name__, h)
                    break
                except Exception:
                    print(f"Exception in mesh generation ({fun.__name__}).")
                    for key, val in vals.items():
                        data[key].append(numpy.nan)
                else:
                    for key, val in vals.items():
                        data[key].append(val)

                data["h"].append(h)
                data["axpy_time"].append(_measure_axpy(vals["num_nodes"]))
                progress.update(task1, advance=1)

            # merge with original file content
            content["data"][fun.__name__] = data

            # store data in files
            with open(module_name.lower() + ".json", "w") as f:
                now = datetime.datetime.utcnow().replace(microsecond=0).isoformat()
                json.dump(
                    {
                        "name": module_name,
                        "date": now,
                        "packages": module.packages,
                        "data": content["data"],
                    },
                    f,
                    indent=2,
                )

            progress.update(task0, advance=1)
    print()


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
            x = numpy.array(data["data"][domain]["num_nodes"])
            y = numpy.array(data["data"][domain]["time"])
            idx = numpy.argsort(x)
            plt.loglog(x[idx], y[idx], color=module.colors[0], label=label)

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
            x = numpy.array(data["data"][domain]["num_nodes"])
            y0 = numpy.array(data["data"][domain]["quality_avg"])
            y1 = numpy.array(data["data"][domain]["quality_min"])
            idx = numpy.argsort(x)
            plt.semilogx(
                x[idx], y0[idx], linestyle="-", color=module.colors[0], label=label
            )
            plt.semilogx(x[idx], y1[idx], linestyle="--", color=module.colors[1])
    dufte.legend()
    plt.xlabel("num points")
    plt.title("cell quality, avg  and min (dashed)")
    plt.ylim(0.0, 1.0)
    plt.savefig(f"{domain}-quality.svg", transparent=True, bbox_inches="tight")
    # plt.show()
    plt.close()

    # num nodes vs. energy
    for module in modules:
        name = module.packages[0][0]
        json_filename = name.lower() + ".json"
        with open(json_filename) as f:
            data = json.load(f)
        label = ", ".join(" ".join(package) for package in module.packages)
        if domain in data["data"]:
            x = numpy.array(data["data"][domain]["num_nodes"])
            idx = numpy.argsort(x)
            try:
                y = numpy.array(data["data"][domain]["energy"])
            except KeyError:
                pass
            else:
                plt.loglog(
                    x[idx], y[idx], linestyle="-", color=module.colors[0], label=label
                )
    dufte.legend()
    plt.xlabel("num points")
    plt.title("energy (lower is better)")
    plt.savefig(f"{domain}-energy.svg", transparent=True, bbox_inches="tight")
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
            x = numpy.array(data["data"][domain]["num_nodes"])
            y = numpy.array(data["data"][domain]["num_poisson_steps"])
            idx = numpy.argsort(x)
            plt.semilogx(x[idx], y[idx], color=module.colors[0], label=label)
    dufte.legend()
    plt.xlabel("num points")
    poisson_tol = 1.0e-10
    plt.title(f"number of CG steps for the Poisson problem (tol={poisson_tol:.1e})")
    plt.ylim(0.0)
    plt.savefig(f"{domain}-poisson.svg", transparent=True, bbox_inches="tight")
    # plt.show()
    plt.close()


def get_poisson_steps_dolfin(pts, cells, tol):
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


def get_poisson_steps_scikitfem(points, cells, tol):
    from skfem import (
        MeshTri,
        MeshTet,
        InteriorBasis,
        ElementTriP1,
        ElementTetP1,
        asm,
        condense,
    )
    from skfem.models.poisson import laplace, unit_load
    import krypy

    if cells.shape[1] == 3:
        mesh = MeshTri(points.T, cells.T)
        e = ElementTriP1()
    else:
        assert cells.shape[1] == 4
        mesh = MeshTet(points.T, cells.T)
        e = ElementTetP1()

    basis = InteriorBasis(mesh, e)

    # assemble
    A = asm(laplace, basis)
    b = asm(unit_load, basis)

    A, b = condense(A, b, I=mesh.interior_nodes(), expand=False)
    _, info = krypy.cg(A, b, tol=tol)

    return info.iter


if __name__ == "__main__":
    for module in modules:
        update_data_file(module)
    for domain, _ in domains_h:
        create_plots(domain)
