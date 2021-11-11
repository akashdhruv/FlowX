import flowx
import time

# Define grid parameters
nx, ny = 448, 448
xmin, xmax = -1, 1
ymin, ymax = -1, 1
xblocks, yblocks = 32,32

# Create flowx objects
varlist = ['phi']
grid = flowx.domain.Grid("cell-centered",varlist,nx,ny,xmin,xmax,ymin,ymax,xblocks,yblocks)
                 
particle_info = dict(input='HDF5', file='sm_tshape.h5')
particle = flowx.domain.Particles(particle_info,xmin,xmax,ymin,ymax,scalars=None)

search_options = dict(monitor=True, nthreads=8, backend='loky')

start = time.time()
ites = flowx.imbound.utils.shapely_search(grid,particle,'phi',search_options)
end = time.time()

print("Search time: ",end-start,"Iterations: ",ites)

start = time.time()
grid.fill_guard_cells('phi')
end = time.time()

print("Fill guardcells time: ",end-start)
