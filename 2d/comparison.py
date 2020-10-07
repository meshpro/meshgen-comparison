import time

import meshplex
import numpy
from rich.progress import Progress

# import cgal_examples
import dmsh_examples
import meshpy_examples
import meshzoo_examples
import pygmsh_examples


def disk():
    H = numpy.logspace(0, -2, num=10)
    times = []
    quality_min = []
    quality_avg = []
    num_cells = []
    modules = {
        # "cgal": cgal_examples,
        "dmsh": dmsh_examples,
        "gmsh": pygmsh_examples,
        "meshpy": meshpy_examples,
        "meshzoo": meshzoo_examples,
    }
    with Progress() as progress:
        task1 = progress.add_task("Overall", total=len(H))
        task2 = progress.add_task("Modules", total=len(modules))
        for h in H:
            times.append([])
            quality_min.append([])
            quality_avg.append([])
            num_cells.append([])
            progress.update(task2, completed=0)
            for module in modules.values():
                tic = time.time()
                points, cells = module.disk(h)
                toc = time.time()

                mesh = meshplex.MeshTri(points, cells)

                times[-1].append(toc - tic)
                quality_min[-1].append(numpy.min(mesh.cell_quality))
                quality_avg[-1].append(numpy.average(mesh.cell_quality))
                num_cells[-1].append(mesh.node_coords.shape[0])
                progress.update(task2, advance=1)
            progress.update(task1, advance=1)

    times = numpy.array(times)
    quality_min = numpy.array(quality_min)
    quality_avg = numpy.array(quality_avg)
    num_cells = numpy.array(num_cells)

    print(times)
    print()
    print(quality_avg)
    print()
    print(quality_min)
    print()
    print(num_cells)



if __name__ == "__main__":
    disk()
