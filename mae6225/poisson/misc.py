import numpy

def get_error(grid,eror,ivar,asol):
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
