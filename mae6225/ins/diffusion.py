"""Routine to compute diffusion terms."""

import numpy


def diffusion(grid, ivar, alpha):
    """Compute the diffusion terms of the variable tagged "ivar".

    Arguments
    ---------
    grid : grid object
        Grid containing data for a given stencil.
    ivar : string
        Name of the grid variable to be operated on.
    alpha : float
        Diffusion coefficient.

    Returns
    -------
    D : numpy.ndarray
        Diffusion terms as an array of floats.

    """
    f = grid.get_values(ivar)

    dx, dy = grid.dx, grid.dy

    D = alpha * ((f[2:, 1:-1] - 2 * f[1:-1, 1:-1] + f[:-2, 1:-1]) / dx**2 +
                 (f[1:-1, 2:] - 2 * f[1:-1, 1:-1] + f[1:-1, :-2]) / dy**2)

    return D
