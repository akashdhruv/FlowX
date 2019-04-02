"""Module with I/O functions."""

import numpy
from matplotlib import pyplot


def plot_contour(grid, var):
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
    pyplot.contourf(X, Y, grid.get_values(var).transpose())
    pyplot.colorbar(label=var)
    pyplot.axis('scaled', adjustable='box')
    pyplot.xlim(X.min(), X.max())
    pyplot.ylim(Y.min(), Y.max())
    pyplot.show()
