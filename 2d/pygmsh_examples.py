import pygmsh


def disk(h):
    with pygmsh.geo.Geometry() as geom:
        geom.add_circle(
            [0.0, 0.0, 0.0],
            1.0,
            mesh_size=h,
            num_sections=4,
            compound=True,
        )
        mesh = geom.generate_mesh()
        mesh.prune_z_0()
    return mesh.points, mesh.get_cells_type("triangle")
