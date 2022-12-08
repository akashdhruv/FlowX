from shapely.geometry import Point, Polygon
from boxkit.library import Action
from boxkit.library import Block
import numpy


def shapely_search(grid, particle, ibmf, options={}):
    """
    Approximate nearest neighbor search

    grid : grid object
         flowx Grid object

    particle : particle object
         flowx Particle object

    ibmf : string
         level set variable

    options : dictionary of options
    """
    points = particle.x[1:, :]
    polygon = Polygon(numpy.ndarray.tolist(points))

    _search_block.monitor = options["monitor"]
    _search_block.nthreads = options["nthreads"]
    _search_block.backend = options["backend"]

    list_count = _search_block(grid.blocklist, polygon, ibmf)
    iter_count = sum(list_count)

    return iter_count


@Action(unit=Block)
def _search_block(unit, polygon, ibmf):
    """
    Block search operation
    """
    iter_count = 0

    for k in range(unit.nzb + 2 * unit.zguard):
        for j in range(unit.nyb + 2 * unit.yguard):
            for i in range(unit.nxb + 2 * unit.xguard):
                unit[ibmf][k, j, i] = (
                    2 * Point([unit.x[i], unit.y[j]]).within(polygon) - 1
                ) * polygon.exterior.distance(Point([unit.x[i], unit.y[j]]))
                iter_count += 1

    return iter_count
