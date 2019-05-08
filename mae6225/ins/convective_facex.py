"""Routine to compute the convective terms at x-faces."""

import numpy


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
    u = gridx.get_values(ivar)
    v = gridy.get_values(ivar)
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
