"""Module with helper functions for the Poisson solver."""

import numpy


def get_error(grid, eror, ivar, asol):
    """Compute the error between the numerical and analytical solutions.

    The error is defined as the absolute difference between the two solutions.

    Arguments
    ---------
    grid : Grid object
        Grid containing the data.
    eror: string
        Name of the grid variable of the error.
    ivar: string
        Name of the grid variable of the numerical solution.
    asol: string
        Name of the grid variable of the analytical solution.

    """
    i_eror, i_ivar, i_asol = grid.get_variable_indices([eror, ivar, asol])
    grid.data[:, :, i_eror] = numpy.abs(grid.data[:, :, i_ivar] -
                                        grid.data[:, :, i_asol])
