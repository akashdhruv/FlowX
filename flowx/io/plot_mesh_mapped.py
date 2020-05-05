"""Module with I/O functions."""

import numpy
from matplotlib import pyplot

def plot_mesh_mapped(grid, particle):
    """Plot the mesh on a meshgrid.

    Arguments
    ---------
    grid : Grid object
        Grid containing the data.
    var : string
        Name of the variable to plot.

    """
    X, Y = numpy.meshgrid((grid.x[:-1]+grid.x[1:])/2, (grid.y[:-1]+grid.y[1:])/2)
    Xc, Yc = numpy.meshgrid(grid.x[1:-1], grid.y[1:-1])

    nodesC = numpy.vstack((particle.x[1:],particle.x[1]))

    pyplot.rc('font', family='serif', size=16)
    pyplot.figure()
    pyplot.xlabel('x')
    pyplot.ylabel('y')
    pyplot.plot(X.T, Y.T, 'g')
    pyplot.plot(X, Y, 'g')
    pyplot.scatter(Xc, Yc)
    pyplot.plot(nodesC[:,0], nodesC[:,1])
    pyplot.plot(nodesC[:,0], nodesC[:,1], marker='.', markevery=5) 
    pyplot.axis('scaled', adjustable='box')
    pyplot.xlim(X.min(), X.max())
    pyplot.ylim(Y.min(), Y.max())
    pyplot.show()
