"""User defined module for simulation"""

import numpy

def get_analytical(grid,asol):
	"""
	Function to calculate exact solution of the variable phi

	Arguments
	---------
	grid : object of class Grid

	asol : scalar integer
	     index where the spatial variable is stored within grid data structure

    
	"""
	
	X,Y = numpy.meshgrid(grid.x_center,grid.y_center)
	X.transpose()
	Y.transpose()

	if(grid.bc_type == 'Neumann' or grid.bc_type == 'neumann'):
		values = numpy.cos(2*numpy.pi*X)*numpy.cos(2*numpy.pi*Y)
	else:
		values = numpy.sin(2*numpy.pi*X)*numpy.sin(2*numpy.pi*Y)

	grid.set_values(asol,values)

	return

def get_RHS(grid,rvar):
	"""
	Function to calculate the RHS of the poisson equation

	Arguments
	---------

	grid : object of class Grid

	rvar : scalar integer
	     index where the spatial variable is stored within the grid data structure

	"""

	X,Y = numpy.meshgrid(grid.x_center,grid.y_center)
	X.transpose()
	Y.transpose()

	if(grid.bc_type == 'Neumann' or grid.bc_type == 'neumann'):
		values = -8*numpy.pi*numpy.pi*numpy.cos(2*numpy.pi*X)*numpy.cos(2*numpy.pi*Y)
	else:
		values = -8*numpy.pi*numpy.pi*numpy.sin(2*numpy.pi*X)*numpy.sin(2*numpy.pi*Y)

	grid.set_values(rvar,values)

	return
