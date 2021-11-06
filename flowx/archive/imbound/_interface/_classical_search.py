import numpy
from numba import jit

def classical_search(x, y, points, nx, ny, np, options):

    phi = numpy.zeros((nx,ny), dtype=float)

    iter_count = _jit_classical_search(x, y, phi, points, nx, ny, np)

    return iter_count, phi

@jit(nopython=True)
def _jit_classical_search(x, y, phi, points, nx, ny, np):
  
    iter_count = 0
   
    for i in range(nx):
        for j in range(ny):

            countit = 0
            dist = 1E13

            for p in range(np):

                # Tag point A
                PA = points[p,:]

                # Tag point B
                if(p == np-1):
                    PB = points[0,:]

                else:
                    PB = points[p+1,:]

                # Tag P1 with grid point co-ordinates
                P1 = numpy.array([x[i,j], y[i,j]])

                # Find intersection of point P1 to node PA-PB
                u = ((P1[0]-PA[0])*(PB[0]-PA[0]) + (P1[1]-PA[1])*(PB[1]-PA[1]))/(((PB[0]-PA[0])**2)+((PB[1]-PA[1])**2))

                # If u is outside the node then change u to point to either of the edges
                if (u < 0.):
                    u = 0.0
                elif (u > 1.):
                    u = 1.0

                # Find the point on the line segment with the shortest distance to P1
                # (If the normal hits the line outside the line segment it is
                # reassigned to hit the closer endpoint.)
                P0 = PA + (PB - PA)*u

                # Determine the quadrant and angle for the "normal"
                # (If to the left or right of the line segment the vector with the 
                #  shortest distance to the line segment will not be perpendicular)
                if (abs(P0[0]-PA[0]) < 1e-13 and abs(P0[1]-PA[1]) < 1e-13):
                    v1 = P1 - P0
                    v2 = P0 - PB
                else:
                    v1 = P1 - P0
                    v2 = PA - P0

                dist = min(dist,numpy.sqrt(v1[0]**2 + v1[1]**2))

                # Find if the horizontal ray on right-side intersects with body
                miny = min(PA[1],PB[1])
                maxy = max(PA[1],PB[1])

                if(y[i,j] > miny and y[i,j] < maxy):

                    # Method #1 use ratios to divide the current panel using
                    # y intersection and find x
                    mratio = PA[1] - y[i,j]
                    nratio = y[i,j] - PB[1]
                    xit = (mratio*PB[0] + nratio*PA[0])/(mratio + nratio)

                    # Method #2 use the equation of line instead
                    #mratio = (PB[1]-PA[1])/(PB[0]-PA[0])          
                    #xit = PA[0] + (y[i,j] - PA[1])/mratio

                    if(xit >= x[i,j]): countit += 1

                iter_count += 1
                
            phi[i,j] = dist
            if(numpy.mod(countit,2) == 0): phi[i,j] = -phi[i,j]

    return iter_count
