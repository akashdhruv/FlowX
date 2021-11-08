"""Module with I/O functions."""

import numpy
from matplotlib import pyplot


def plot_contour(grid, ivar, levels=None):
    """Plot the filled contour of a ivariable on a meshgrid.

    Arguments
    ---------
    grid : Grid object
        Grid containing the data.
    ivar : string
        Name of the ivariable to plot.

    """
    X, Y = numpy.meshgrid(grid.x, grid.y)

    pyplot.rc('font', family='serif', size=16)
    pyplot.figure()
    pyplot.xlabel('x')
    pyplot.ylabel('y')

    if not levels:
        pyplot.contourf(X, Y, grid[ivar])
    else:
        pyplot.contourf(X, Y, grid[ivar], levels)
   
    pyplot.colorbar(label=ivar)
    pyplot.axis('scaled') #adjustable='box')
    pyplot.xlim(X.min(), X.max())
    pyplot.ylim(Y.min(), Y.max())
    #pyplot.show()
