import numpy as np
from lattice_gauge_theory import lattice, lattice_site
import re

def square_lattice(a, b, spacing=1.0):
    """
    Generate a square lattice.

    Args:
        a (int): Number of lattice sites along x.
        b (int): Number of lattice sites along y.
        spacing (float): Distance between lattice sites. (Default: 1.0)

    Returns:
        (Lattice): The initialized lattice.

    Notes:
        The returned lattice is 3D periodic, but all sites and edges lie in the
        xy-plane.
    """
    grid = np.array(list(range(1, a*b + 1))).reshape(a, b, order='F')
    it = np.nditer(grid, flags=['multi_index'])
    sites = []
    while not it.finished:
        x, y = it.multi_index
        r = np.array([x*spacing, y*spacing, 0.0])
        neighbors = [
            np.roll(grid, +1, axis=0)[x,y],
            np.roll(grid, -1, axis=0)[x,y],
            np.roll(grid, +1, axis=1)[x,y],
            np.roll(grid, -1, axis=1)[x,y],
        ]
        sites.append(lattice_site.Site(int(it[0]), r, neighbors, 0.0, 'L'))
        it.iternext()
    return lattice.Lattice(sites,
                           cell_lengths=np.array([a, b, 0.0])*spacing,
                           dim=2)

def cubic_lattice(a, b, c, spacing):
    """
    Generate a cubic lattice.

    Args:
        a (int): Number of lattice sites along x.
        b (int): Number of lattice sites along y.
        c (int): Number of lattice sites along z.
        spacing (float): Distance between lattice sites. (Default: 1.0)

    Returns:
        (Lattice): The initialized lattice.
    """
    grid = np.array(list(range(1, a*b*c + 1))).reshape(a, b, c, order='F')
    it = np.nditer(grid, flags=['multi_index'])
    while not it.finished:
        x, y, z = it.multi_index
        r = np.array([x, y, z]) * spacing
        neighbors = [
            np.roll(grid, +1, axis=0)[x,y,z],
            np.roll(grid, -1, axis=0)[x,y,z],
            np.roll(grid, +1, axis=1)[x,y,z],
            np.roll(grid, -1, axis=1)[x,y,z],
            np.roll(grid, +1, axis=2)[x,y,z],
            np.roll(grid, -1, axis=2)[x,y,z],
        ]
        sites.append(lattice_site.Site(int(it[0]), r, neighbors, 0.0, 'L'))
        it.iternext()
    return lattice.Lattice(sites,
                           cell_lengths=np.array([a, b, c]) *spacing,
                           dim=3)


