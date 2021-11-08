"""Module with I/O functions."""

import numpy
from matplotlib import pyplot


def plot_vector(gridx, gridy, var):
    """Plot the vector for given x, y data.

    Arguments
    ---------
    gridx : Grid object
        Grid containing the data on x-face.
    gridy : Grid object
        Grid containing the data on y-face.
    var : string
        Name of the variable to plot.

    """
    X, Y = numpy.meshgrid(gridx.x, gridy.y)

    u = gridx.get_values(var)
    v = gridy.get_values(var)

    U = (u[:, 1:] + u[:, :-1]) / 2
    V = (v[1:, :] + v[-1:, :]) / 2

    pyplot.rc('font', family='serif', size=16)
    pyplot.figure()
    pyplot.xlabel('x')
    pyplot.ylabel('y')
    pyplot.quiver(X, Y, U.transpose(), V.transpose())
    pyplot.axis('scaled')#, adjustable='box')
    pyplot.xlim(X.min(), X.max())
    pyplot.ylim(Y.min(), Y.max())
    pyplot.show()
