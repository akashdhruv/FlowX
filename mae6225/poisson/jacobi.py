import numpy

def solve_jacobi(grid,ivar,rvar,max_iterations=3000,tol=1e-9):
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
	phi_old           = numpy.copy(grid.data[:,:,grid.center_vars[ivar]])

	while iteration_counter <= max_iterations:

		grid.data[1:-1,1:-1,grid.center_vars[ivar]] = ((phi_old[1:-1,2:]/(grid.dy**2))  + \
		                                               (phi_old[1:-1,:-2]/(grid.dy**2)) + \
					                       (phi_old[2:,1:-1]/(grid.dx**2))  + \
				                               (phi_old[:-2,1:-1]/(grid.dx**2)) - \
				                                grid.data[1:-1,1:-1,grid.center_vars[rvar]])       * \
                                                               (1/((1/(grid.dx**2))             + \
                                                               (1/(grid.dy**2))             + \
                                                               (1/(grid.dx**2))             + \
                                                               (1/(grid.dy**2))))

		grid.fill_guard_cells([ivar])

		residual = numpy.sqrt(numpy.sum((grid.data[:,:,grid.center_vars[ivar]]-phi_old)**2)/((grid.nx+2)*(grid.ny+2)))

		if(residual < tol and residual != 0.0):
			break			

		iteration_counter = iteration_counter + 1
		phi_old           = numpy.copy(grid.data[:,:,grid.center_vars[ivar]])

	return iteration_counter, residual
