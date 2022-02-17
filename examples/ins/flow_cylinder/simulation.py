"""User-defined module for simulation."""

import numpy
import flowx


def set_initial_velocity(gridc, gridx, gridy, ivar, pres):
    """Set the initial velocity field.

    The x- and y-components of the velocity are set to 1.0 and 0.0,
    respectively.

    Arguments
    ---------
    gridx : flowx.Grid object
        Grid containing x-face data.
    gridy : flowx.Grid object
        Grid containing y-face data.
    ivar : string
        Name of the velocity variable on the grid.

    """

    u = gridx[ivar]
    v = gridy[ivar]
    p = gridc[pres]

    u[:, :] = 1.0
    v[:, :] = 0.0
    p[:, :] = 0.0

    return
