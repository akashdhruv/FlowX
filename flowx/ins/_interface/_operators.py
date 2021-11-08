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
    f = grid[ivar]

    dx, dy = grid.dx, grid.dy

    D = alpha * ((f[2:, 1:-1] - 2 * f[1:-1, 1:-1] + f[:-2, 1:-1]) / dx**2 +
                 (f[1:-1, 2:] - 2 * f[1:-1, 1:-1] + f[1:-1, :-2]) / dy**2)

    return D

def convective_facex(gridx, gridy, ivar):
    """Convection operator for the x-face grid.

    Arguments
    ---------
    gridx : grid object (x-direction)
        Grid containing data in x-direction.
    gridy : grid object (y-direction)
        Grid containing data in y-direction.
    ivar: string
        Name of the grid variable to be operated on.

    Returns
    -------
    F : numpy.ndarray
        Convective terms in the x direction as an array of floats.

    """
    u = gridx[ivar]
    v = gridy[ivar]
    dx, dy = gridx.dx, gridy.dy

    u_P = u[1:-1, 1:-1]
    u_W = u[:-2, 1:-1]
    u_E = u[2:, 1:-1]
    u_S = u[1:-1, :-2]
    u_N = u[1:-1, 2:]

    v_sw = v[1:-2, :-1]
    v_se = v[2:-1, :-1]
    v_nw = v[1:-2, 1:]
    v_ne = v[2:-1, 1:]

    F = - (((u_P + u_E)**2 - (u_W + u_P)**2) / (4 * dx) +
           ((u_P + u_N) * (v_nw + v_ne) -
            (u_S + u_P) * (v_sw + v_se)) / (4 * dy))

    return F

def convective_facey(gridx, gridy, ivar):
    """Convection operator for the the x-face grid.

    Arguments
    ---------
    gridx : grid object (x-direction)
        Grid containing data in x-direction.
    gridy : grid object (y-direction)
        Grid containing data in y-direction.
    ivar : string
        Name of the grid variable to be operated on.

    Returns
    -------
    F : numpy.ndarray
        Convective terms in the y direction as an array of floats.

    """
    u = gridx[ivar]
    v = gridy[ivar]
    dx, dy = gridx.dx, gridy.dy

    v_P = v[1:-1, 1:-1]
    v_W = v[:-2, 1:-1]
    v_E = v[2:, 1:-1]
    v_S = v[1:-1, :-2]
    v_N = v[1:-1, 2:]

    u_sw = u[:-1, 1:-2]
    u_se = u[1:, 1:-2]
    u_nw = u[:-1, 2:-1]
    u_ne = u[1:, 2:-1]

    F = -(((u_se + u_ne) * (v_P + v_E) -
           (u_sw + u_nw) * (v_W + v_P)) / (4 * dx) +
          ((v_P + v_N)**2 - (v_S + v_P)**2) / (4 * dy))

    return F
