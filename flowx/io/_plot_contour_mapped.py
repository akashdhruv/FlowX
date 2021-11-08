"""Module with I/O functions."""

import numpy
from matplotlib import pyplot


def plot_contour_mapped(grid, ivar, particle, levels=None):
    """Plot the filled contour of a variable on a meshgrid.

    Arguments
    ---------
    grid : Grid object
        Grid containing the data.
    ivar : string
        Name of the variable to plot.

    """
    X, Y = numpy.meshgrid(grid.x, grid.y)

    nodesC = numpy.vstack((particle.x[1:],particle.x[1]))

    pyplot.rc('font', family='serif', size=16)
    pyplot.figure()
    pyplot.xlabel('x')
    pyplot.ylabel('y')

    if not levels:
        pyplot.contourf(X, Y, grid[ivar])
    else:
        pyplot.contourf(X, Y, grid[ivar], levels)

    pyplot.plot(nodesC[:,0], nodesC[:,1])
    pyplot.plot(nodesC[:,0], nodesC[:,1], marker='.', markevery=5)
   
    pyplot.colorbar(label=ivar)
    pyplot.axis('scaled')#, adjustable='box')
    pyplot.xlim(X.min(), X.max())
    pyplot.ylim(Y.min(), Y.max())
    pyplot.show()
