"""User defined module for simulation."""

import numpy


def get_analytical(grid, asol):
    """Compute and set the analytical solution.

    Arguments
    ---------
    grid : mae6225.Grid object
        Grid containing data.
    asol : string
        Name of the variable on the grid.

    """
    X, Y = numpy.meshgrid(grid.x_center, grid.y_center)
    values = numpy.sin(2 * numpy.pi * X) * numpy.sin(2 * numpy.pi * Y)
    grid.set_values(asol, values.transpose())


def get_rhs(grid, rvar):
    """Compute and set the right-hand side of the Poisson system.

    Arguments
    ---------
    grid : mae6225.Grid object
        Grid containing data.
    rvar : string
        Name of the variable on the grid.

    """
    X, Y = numpy.meshgrid(grid.x_center, grid.y_center)
    values = (-8 * numpy.pi**2 *
              numpy.sin(2 * numpy.pi * X) * numpy.sin(2 * numpy.pi * Y))
    grid.set_values(rvar, values.transpose())
