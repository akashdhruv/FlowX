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
    fx = gridx.get_values(ivar)
    fy = gridy.get_values(ivar)

    dx = gridx.dx
    dy = gridy.dy

    val_xE = fx[2:, 1:-1]
    val_xW = fx[:-2, 1:-1]
    val_xP = fx[1:-1, 1:-1]
    val_xN = fx[1:-1, 2:]
    val_xS = fx[1:-1, :-2]

    val_yNW = fy[1:-2, 1:]
    val_yNE = fy[2:-1, 1:]
    val_ySW = fy[1:-2, :-1]
    val_ySE = fy[2:-1, :-1]

    F = (-((val_xP + val_xE)**2 - (val_xP + val_xW)**2) / (4 * dx) -
         (((val_xP + val_xN) * (val_yNE + val_yNW)) / (4 * dy) -
          ((val_xP + val_xS) * (val_ySE + val_ySW)) / (4 * dy)))

    return F
