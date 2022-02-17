"""Module with I/O functions."""

import numpy
import matplotlib.pyplot as pyplot
import matplotlib.ticker as ticker


def plot_contour(grid, ivar, auxvars=[], particle=None, scalars=None, levels=None):
    """Plot the filled contour of a variable on a meshgrid.

    Arguments
    ---------
    grid : Grid object
        Grid containing the data.
    ivar : string
        Name of the variable to plot.
    """
    pyplot.rc("font", family="serif", size=16)
    pyplot.figure()
    pyplot.xlabel("x")
    pyplot.ylabel("y")

    ivarmin = numpy.min(grid[ivar])
    ivarmax = numpy.max(grid[ivar])

    if not levels:
        levels = numpy.linspace(ivarmin, ivarmax, 10)

    for block in grid.blocklist:
        xmesh, ymesh = numpy.meshgrid(block.x, block.y)

        if len(levels) > 1:
            pyplot.contourf(xmesh, ymesh, block[ivar][0, :, :], levels=levels)
        else:
            pyplot.contour(
                xmesh,
                ymesh,
                block[ivar][0, :, :],
                levels=levels,
                linestyles="solid",
                linewidths=1.0,
                colors="black",
            )

        for var in auxvars:
            pyplot.contour(
                xmesh,
                ymesh,
                block[var][0, :, :],
                colors="blue",
                linestyles="solid",
                linewidths=1.0,
            )

    if particle:
        nodes = numpy.vstack((particle.x[1:], particle.x[1]))
        pyplot.plot(nodes[:, 0], nodes[:, 1], color="maroon")
        pyplot.plot(
            nodes[:, 0],
            nodes[:, 1],
            marker=".",
            markersize=8,
            markevery=4,
            color="maroon",
        )

    if len(levels) > 1:
        pyplot.colorbar(label=ivar, format=ticker.FuncFormatter(_fmt))
    pyplot.axis("scaled")
    pyplot.xlim(grid.xmin, grid.xmax)
    pyplot.ylim(grid.ymin, grid.ymax)
    # pyplot.savefig('./images/grid%d.png' % int(scalars.nstep/100))
    # pyplot.show()


def _fmt(x, pos):
    a, b = "{:.0e}".format(x).split("e")
    b = int(b)
    return r"${} \times 10^{{{}}}$".format(a, b)
