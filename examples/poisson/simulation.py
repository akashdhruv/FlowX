"""User defined module for simulation."""

import numpy


def get_analytical(grid, asol, user_bc):
    """Compute and set the analytical solution.

    Arguments
    ---------
    grid : flowx.Grid object
        Grid containing data.
    asol : string
        Name of the variable on the grid.

    """
    for block in grid.blocklist:
        xmesh, ymesh = numpy.meshgrid(block.x, block.y)

        if(user_bc == 'dirichlet'):
            values = numpy.sin(2 * numpy.pi * xmesh) * numpy.sin(2 * numpy.pi * ymesh)
        else:
            values = numpy.cos(2 * numpy.pi * xmesh) * numpy.cos(2 * numpy.pi * ymesh)

        block[asol] = values

def get_rhs(grid, rvar, user_bc):
    """Compute and set the right-hand side of the Poisson system.

    Arguments
    ---------
    grid : flowx.Grid object
        Grid containing data.
    rvar : string
        Name of the variable on the grid.

    """
    for block in grid.blocklist:
        xmesh, ymesh = numpy.meshgrid(block.x, block.y)

        if(user_bc == 'dirichlet'):
            values = (-8 * numpy.pi**2 *
                      numpy.sin(2 * numpy.pi * xmesh) * numpy.sin(2 * numpy.pi * ymesh))
        else:
            values = (-8 * numpy.pi**2 *
                      numpy.cos(2 * numpy.pi * xmesh) * numpy.cos(2 * numpy.pi * ymesh))

        block[rvar] = values
