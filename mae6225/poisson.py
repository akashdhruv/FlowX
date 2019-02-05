"Module for Poisson solver unit"

import numpy

def getAnalytical(grid,asol):
	"""
	Function to calculate exact solution of the variable phi

	Arguments
	---------
	grid : object of class Grid

	asol : scalar integer
	     index where the spatial variable is stored within grid data structure

    
	"""
	
	for i in range(grid.nx+2):
		for j in range(grid.ny+2):
			value = numpy.cos(2*numpy.pi*grid.x_center[i])*numpy.cos(2*numpy.pi*grid.y_center[j])
			grid.set_value(asol,i,j,value)
	return

def getRHS(grid,rvar):
	"""
	Function to calculate the RHS of the poisson equation

	Arguments
	---------

	grid : object of class Grid

	rvar : scalar integer
	     index where the spatial variable is stored within the grid data structure

	"""

	for i in range(grid.nx+2):
		for j in range(grid.ny+2):
			value = -8*numpy.pi*numpy.pi*numpy.cos(2*numpy.pi*grid.x_center[i])*numpy.cos(2*numpy.pi*grid.y_center[j])
			grid.set_value(rvar,i,j,value)
	return

def solveJacobi(grid,ivar,rvar,max_iterations,tol):
	"""
	Function to solve the Poisson equation using iterative (Jacobi) method

	Arguments
	---------

	grid : object of class Grid

	ivar : scalar integer
	     index to store numerical solution

	rvar : scalar integer
	     index where the right hand side is stored

	max_iterations : scalar integer
		       limit on number of iterations

	tol : scalar float
	    tolerance for convergence check

	-----------------
	Returned variables
	----------------

	iteration_counter : scalar integer
			  total number of iterations until convergence

	residual : scalar float
		 residual error between after completion

	---------------
	Local variables
	---------------

	phi_old : float array
		variable to store numerical solution from the previous iteration

	"""

	iteration_counter = 1
	phi_old           = numpy.copy(grid.data[:,:,ivar])

	while iteration_counter <= max_iterations:

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

		residual = numpy.sqrt(numpy.sum((grid.data[:,:,ivar]-phi_old)**2)/((grid.nx+2)*(grid.ny+2)))

		if(residual < tol and residual != 0.0):
			break			

		iteration_counter = iteration_counter + 1
		phi_old           = numpy.copy(grid.data[:,:,ivar])

	return iteration_counter, residual

def getError(grid,eror,ivar,asol):
	"""
	Function to calculate absoulte error between numerical and exact solution

	Arguments
	---------

	grid : object of class Grid

	eror : scalar integer
	     index to store the absolute error within the grid data structure

	ivar : scalar integer
	     index where the numerical solution is stored

	asol : scalar integer
	     index where the exact solution is stored

	"""

	grid.data[:,:,eror] = numpy.absolute(grid.data[:,:,ivar]-grid.data[:,:,asol])

	return

def applyBC(phi):	

	"""
	Function to apply boundary conditions

	Arguments
	---------

	phi : float array
	    variable on which the BC needs to be applied
	"""

	# Homogeneous Neumann BC
	phi[:,0]  = phi[:,1]
	phi[:,-1] = phi[:,-2]
	phi[0,:]  = phi[1,:]
	phi[-1,:] = phi[-2,:]

	# Homogeneous Dirichlet BC
	#phi[:,0]  = -phi[:,1]
	#phi[:,-1] = -phi[:,-2]
	#phi[0,:]  = -phi[1,:]
	#phi[-1,:] = -phi[-2,:]

	return
