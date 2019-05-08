"""Routine to compute the convective terms at y-faces."""

import numpy


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
    u = gridx.get_values(ivar)
    v = gridy.get_values(ivar)
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
