import numpy
from numba import jit
from shapely.geometry import Point, Polygon

def shapely_search(x, y, points, nx, ny, np, options):

    iter_count = 0
 
    phi = numpy.zeros((ny,nx),dtype=float) 
    max_panel_length = options['max_panel_length']
    polygon = Polygon(numpy.ndarray.tolist(points))

    for j in range(ny):
        for i in range(nx):    
            phi[j,i] = (2*Point([x[j,i], y[j,i]]).within(polygon) - 1)*polygon.exterior.distance(Point([x[j,i], y[j,i]]))    
            iter_count += 1
    
    return iter_count, phi
