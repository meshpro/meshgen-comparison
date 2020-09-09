import pygmsh


def disk(h):
    geom = pygmsh.built_in.Geometry()
    geom.add_circle(
        [0.0, 0.0, 0.0],
        1.0,
        lcar=h,
        num_sections=4,
        compound=True,
    )
    mesh = pygmsh.generate_mesh(geom, prune_z_0=True, verbose=False)
    return mesh.points, mesh.get_cells_type("triangle")
