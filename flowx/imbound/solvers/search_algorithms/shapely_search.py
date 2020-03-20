import numpy
from numba import jit
from shapely.geometry import Point, Polygon

def shapely_search(x, y, phi, points, nx, ny, np, max_panel_length, options):

    iter_count = 0

    polygon = Polygon(numpy.ndarray.tolist(points))

    for i in range(nx):
        for j in range(ny):

            phi[i,j] = (2*Point([x[i,j], y[i,j]]).within(polygon) - 1)*polygon.exterior.distance(Point([x[i,j], y[i,j]]))

            iter_count += 1
 
    return iter_count
