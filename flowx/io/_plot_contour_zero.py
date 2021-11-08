"""Module with I/O functions."""

import numpy
from matplotlib import pyplot


def plot_contour_zero(grid, scalars, ivar, xvar, yvar):
    """Plot the filled contour of a variable on a meshgrid.

    Arguments
    ---------
    grid : Grid object
        Grid containing the data.
    var : string
        Name of the variable to plot.

    """
    X, Y = numpy.meshgrid(grid.x, grid.y)

    pyplot.rc('font', family='serif', size=16)
    pyplot.figure()
    pyplot.xlabel('x')
    pyplot.ylabel('y')

    pyplot.contour(X, Y, grid[xvar], colors='black', linestyles='solid', linewidths=1.0)
    pyplot.contour(X, Y, grid[yvar], colors='black', linestyles='solid', linewidths=1.0)
    pyplot.contour(X, Y, grid[ivar], levels=[0], colors='red', linestyles='solid', linewidths=1.5)
   
    pyplot.axis('scaled') #adjustable='box')
    pyplot.xlim(X.min(), X.max())
    pyplot.ylim(Y.min(), Y.max())
    #pyplot.savefig('./images/grid%d.png' % int(scalars.nstep/100))
    #pyplot.show()
