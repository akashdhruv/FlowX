import numpy
from numba import jit
from shapely.geometry import Point, Polygon

def shapely_search(x, y, phi, nx, ny, particle, options):

    iter_count = 0
 
    max_panel_length = particle.max_panel_length

    points = particle.x[1:,:]

    polygon = Polygon(numpy.ndarray.tolist(points))

    for i in range(nx):
        for j in range(ny):
    
            phi[i,j] = (2*Point([x[i,j], y[i,j]]).within(polygon) - 1)*polygon.exterior.distance(Point([x[i,j], y[i,j]]))
    
            iter_count += 1
    

    return iter_count
