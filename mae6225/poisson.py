import numpy

def getAnalytical(grid,asol):
	for i in range(grid.nx+2):
		for j in range(grid.ny+2):
			value = numpy.cos(2*numpy.pi*grid.x_center[i])*numpy.cos(2*numpy.pi*grid.y_center[j])
			grid.set_value(asol,i,j,value)
	return

def getRHS(grid,rvar):
	for i in range(grid.nx+2):
		for j in range(grid.ny+2):
			value = -8*numpy.pi*numpy.pi*numpy.cos(2*numpy.pi*grid.x_center[i])*numpy.cos(2*numpy.pi*grid.y_center[j])
			grid.set_value(rvar,i,j,value)
	return

def solveJacobi(grid,ivar,rvar):

	max_iteration     = 5000
	iteration_counter = 1
	phi_old           = numpy.copy(grid.data[:,:,ivar])

	while iteration_counter < max_iteration:

		grid.data[1:-1,1:-1,ivar] = ((phi_old[1:-1,2:]/(grid.dy**2))  + \
		                             (phi_old[1:-1,:-2]/(grid.dy**2)) + \
					     (phi_old[2:,1:-1]/(grid.dx**2))  + \
				             (phi_old[:-2,1:-1]/(grid.dx**2)) - \
				             grid.data[1:-1,1:-1,rvar])       * \
                                             (1/((1/(grid.dx**2))             + \
                                                 (1/(grid.dy**2))             + \
                                                 (1/(grid.dx**2))             + \
                                                 (1/(grid.dy**2))))

		applyBC(grid.data[:,:,ivar])

		residual = numpy.sum((grid.data[:,:,ivar]-phi_old)**2)/((grid.nx+2)*(grid.ny+2))
	
		if(residual < 0.000001 and residual != 0.0):
			break

		iteration_counter = iteration_counter + 1
		phi_old           = numpy.copy(grid.data[:,:,ivar])

	return iteration_counter, residual

def getError(grid,eror,ivar,asol): 

	return

def applyBC(phi):	

	phi[:,0]  = phi[:,1]
	phi[:,-1] = phi[:,-2]
	phi[0,:]  = phi[1,:]
	phi[-1,:] = phi[-2,:]

	return
