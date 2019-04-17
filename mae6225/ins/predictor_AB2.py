"""Routine to compute the intermediate velocity."""

import numpy

from .convective_facex import *
from .convective_facey import *
from .diffusion import *


def predictor_AB2(gridx, gridy, ivar, hvar, Re, ifac):
    """Velocity prediction step in x and y direction.

    Arguments
    ---------
    gridx : grid object (x-direction)
        Grid containing data in x-direction.
    gridy : grid object (y-direction)
        Grid containing data in y-direction.
    ivar : string
        Name of the grid variable of the velocity solution.
    hvar : string
        Name of the grid variable to store convective + diffusion terms.
    Re : float
        Reynolds number.
    ifac : float
        Time-step size.

    """
    hx_old = gridx.get_values(hvar)
    hy_old = gridy.get_values(hvar)

    hx_new = (convective_facex(gridx, gridy, ivar) +
                      diffusion(gridx, ivar, 1 / Re)) 
    hy_new = (convective_facey(gridx, gridy, ivar) +
                      diffusion(gridy, ivar, 1 / Re)) 

    u = gridx.get_values(ivar)
    v = gridy.get_values(ivar)

    u[1:-1, 1:-1] = u[1:-1, 1:-1] + ifac * (1.5 * hx_new - 0.5 * hx_old[1:-1, 1:-1])
    v[1:-1, 1:-1] = v[1:-1, 1:-1] + ifac * (1.5 * hy_new - 0.5 * hy_old[1:-1, 1:-1])
    
    hx_old[1:-1,1:-1] = hx_new
    hy_old[1:-1,1:-1] = hy_new

    gridx.fill_guard_cells(ivar)
    gridy.fill_guard_cells(ivar)

    return
