import numpy
from numba import jit
from annoy import AnnoyIndex

def ann_search(x, y, phi, points, nx, ny, np, tracing_tol, options):

    iter_count = 0

    nodesA = points
    nodesB = numpy.vstack((points[1:],points[0]))
    nodesC = (nodesA + nodesB)/2.0
    nodesY = nodesC[:,1]

    tree_size = 2
    trace_size = 1

    tree = AnnoyIndex(tree_size, 'euclidean')
    trace = AnnoyIndex(trace_size, 'euclidean')

    ntrees = options['ntrees']
    nquery_trees = options['nquery_trees']
    nquery_trace = options['nquery_trace']

    for i in range(np): 
        tree.add_item(i,nodesC[i])
        trace.add_item(i,[nodesY[i]])

    tree.build(ntrees)
    trace.build(ntrees)

    for i in range(nx):
        for j in range(ny):

            tree_index = tree.get_nns_by_vector([x[i,j],y[i,j]], nquery_trees)
            dist = 1E13

            trace_index, trace_distance = trace.get_nns_by_vector([y[i,j]], nquery_trace, include_distances=True) 
            countit = 0

            niter_trace = numpy.size(numpy.where(numpy.array(trace_distance) < tracing_tol))

            for p in range(nquery_trees):

                PA = nodesA[tree_index[p]]
                PB = nodesB[tree_index[p]]

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
           
                iter_count += 1

            phi[i,j] = -dist

            for p in range(niter_trace):

                PA = nodesA[trace_index[p]]
                PB = nodesB[trace_index[p]]

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

            if(numpy.mod(countit,2) == 1): phi[i,j] = -phi[i,j]

    return iter_count
