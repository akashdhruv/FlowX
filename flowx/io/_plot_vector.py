"""Module with I/O functions."""

import numpy
from matplotlib import pyplot


def plot_vector(gridx, gridy, ivar):
    """Plot the vector for given x, y data.

    Arguments
    ---------
    gridx : Grid object
        Grid containing the data on x-face.
    gridy : Grid object
        Grid containing the data on y-face.
    ivar : string
        Name of the ivariable to plot.

    """
    X, Y = numpy.meshgrid(gridx.x, gridy.y)

    u = gridx[ivar]
    v = gridy[ivar]

    U = (u[1:, :] + u[:-1, :]) / 2
    V = (v[:, 1:] + v[:, -1:]) / 2

    pyplot.rc('font', family='serif', size=16)
    pyplot.figure()
    pyplot.xlabel('x')
    pyplot.ylabel('y')
    pyplot.quiver(X, Y, U, V)
    pyplot.axis('scaled')#, adjustable='box')
    pyplot.xlim(X.min(), X.max())
    pyplot.ylim(Y.min(), Y.max())
    pyplot.show()
