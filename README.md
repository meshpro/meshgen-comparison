# Mesh Generator Comparison

This repository creates meshes of the same domains with multiple mesh generators and
compares the results. The used mesh generators are

  * [Gmsh](https://gmsh.info/) (via [pygmsh](https://github.com/nschloe/pygmsh)),
  * [CGAL](https://www.cgal.org/) (via [pygalmesh](https://github.com/nschloe/pygalmesh)),
  * [MeshPy](https://github.com/inducer/meshpy),
  * [dmsh](https://github.com/nschloe/dmsh),
  * [meshzoo](https://github.com/nschloe/meshzoo), and
  * [SeismicMesh](https://github.com/krober10nd/SeismicMesh).

For each domain and mesh generator (if applicable), the repostory shows the following data:
   * the mesh generation time
   * the average and minimum cell quality
   * the nummber of CG iterations for the FEM-Poisson-problem on the mesh.

### Disk

blank | <img src="https://nschloe.github.io/meshgen-comparison/disk-times.svg" width="100%"> |
:-------------:|:-----------------:|
Example mesh   |  Generation time  |
| <img src="https://nschloe.github.io/meshgen-comparison/disk-quality.svg" width="100%"> |<img src="https://nschloe.github.io/meshgen-comparison/disk-poisson.svg" width="100%"> |
Cell quality   |  Poisson CG steps  |


### L-shape

### Rectangle with refinement

### Ball

### L-shape in 3D

### Box with refinement
