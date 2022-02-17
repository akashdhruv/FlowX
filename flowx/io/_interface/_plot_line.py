"""Module with I/O functions."""

import numpy
from matplotlib import pyplot


def plot_line(particle):
    """Plot the vector for given x, y data.

    Arguments
    ---------
    """
    nodes = numpy.vstack((particle.x[1:], particle.x[1]))

    pyplot.rc("font", family="serif", size=16)
    pyplot.figure()
    pyplot.xlabel("x")
    pyplot.ylabel("y")
    pyplot.plot(nodes[:, 0], nodes[:, 1])
    pyplot.plot(nodes[:, 0], nodes[:, 1], marker=".", markevery=5)
    pyplot.axis("scaled")
    pyplot.xlim(nodes[:, 0].min(), nodes[:, 0].max())
    pyplot.ylim(nodes[:, 1].min(), nodes[:, 1].max())
    pyplot.show()
