"""Module with I/O functions."""

import numpy
from matplotlib import pyplot

def plot_mesh(grid, particle=None):
    """Plot the mesh on a meshgrid.

    Arguments
    ---------
    grid : Grid object
        Grid containing the data.
    var : string
        Name of the variable to plot.

    """
    xmesh, ymesh   = numpy.meshgrid((grid.x[:-1]+grid.x[1:])/2, (grid.y[:-1]+grid.y[1:])/2)
    xcenter, ycenter = numpy.meshgrid(grid.x[1:-1], grid.y[1:-1])

    pyplot.rc('font', family='serif', size=16)
    pyplot.figure()
    pyplot.xlabel('x')
    pyplot.ylabel('y')
    pyplot.plot(xmesh.T, ymesh.T, 'g')
    pyplot.plot(xmesh, ymesh, 'g')
    pyplot.scatter(xcenter, ycenter)

    if particle:
        nodes = numpy.vstack((particle.x[1:],particle.x[1]))
        pyplot.plot(nodes[:,0], nodes[:,1], color='maroon')
        pyplot.plot(nodes[:,0], nodes[:,1], marker='.', markersize=8, markevery=4, color='maroon')  

    pyplot.axis('scaled')
    pyplot.xlim(grid.xmin, grid.xmax)
    pyplot.ylim(grid.ymin, grid.ymax)
    pyplot.show()
