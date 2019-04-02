"""Routine to compute the final divergence-free velocity."""

import numpy


def corrector(gridc, gridx, gridy, ivar, pvar, ifac):
    """Velocity correction in x and y direction.

    Arguments
    ---------
    gridc : grid object (cell center)
        Grid contaning data in cell center.
    gridx : grid object (x-direction)
        Grid containing data in x-direction.
    gridy : grid object  (y-direction)
        Grid containing data in y-direction.
    ivar : string
        Name of the grid variable of the velocity solution.
    pvar : string
        Name of the grid variable of the pressure solution.
    ifac : float
        Time-step size.

    """
    u = gridx.get_values(ivar)
    v = gridy.get_values(ivar)
    p = gridc.get_values(pvar)

    dx, dy = gridc.dx, gridc.dy

    u[1:-1, 1:-1] = u[1:-1, 1:-1] - ifac * (p[2:-1, 1:-1] - p[1:-2, 1:-1]) / dx
    v[1:-1, 1:-1] = v[1:-1, 1:-1] - ifac * (p[1:-1, 2:-1] - p[1:-1, 1:-2]) / dy

    gridx.fill_guard_cells(ivar)
    gridy.fill_guard_cells(ivar)

    return
