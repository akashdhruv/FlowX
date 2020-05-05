"""Module with I/O functions."""
  
import numpy
from matplotlib import pyplot


def plot_scatter(particles):
    """Plot the vector for given x, y data.

    Arguments
    ---------
    """

    nodesA = numpy.vstack((particles[0].x[1:],particles[0].x[1]))
    nodesB = particles[1].x    

    pyplot.rc('font', family='serif', size=16)
    pyplot.figure()
    pyplot.xlabel('x')
    pyplot.ylabel('y')
    pyplot.plot(nodesA[:,0], nodesA[:,1])
    pyplot.plot(nodesA[:,0], nodesA[:,1], marker='.', markevery=5)
    pyplot.scatter(nodesB[:,0], nodesB[:,1])
    pyplot.axis('scaled', adjustable='box')
    pyplot.axis('equal')
    pyplot.show()
