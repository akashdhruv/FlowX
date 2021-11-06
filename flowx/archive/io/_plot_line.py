"""Module with I/O functions."""
  
import numpy
from matplotlib import pyplot


def plot_line(particle):
    """Plot the vector for given x, y data.

    Arguments
    ---------
    """


    nodesC = numpy.vstack((particle.x[1:],particle.x[1]))

    pyplot.rc('font', family='serif', size=16)
    pyplot.figure()
    pyplot.xlabel('x')
    pyplot.ylabel('y')
    pyplot.plot(nodesC[:,0], nodesC[:,1])
    pyplot.plot(nodesC[:,0], nodesC[:,1], marker='.', markevery=5)
    pyplot.axis('scaled')#, adjustable='box')
    pyplot.xlim(nodesC[:,0].min(), nodesC[:,0].max())
    pyplot.ylim(nodesC[:,1].min(), nodesC[:,1].max())
    pyplot.show()
