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
    xmesh, ymesh = numpy.meshgrid(gridx.x, gridy.y)

    uface = gridx[ivar][0, 0, :, :]
    vface = gridy[ivar][0, 0, :, :]

    umesh = (uface[1:, :] + uface[:-1, :]) / 2
    vmesh = (vface[:, 1:] + vface[:, -1:]) / 2

    pyplot.rc("font", family="serif", size=16)
    pyplot.figure()
    pyplot.xlabel("x")
    pyplot.ylabel("y")
    pyplot.quiver(xmesh, ymesh, umesh, vmesh)
    pyplot.axis("scaled")
    pyplot.xlim(gridx.xmin, gridx.xmax)
    pyplot.ylim(gridy.ymin, gridy.ymax)
    pyplot.show()
