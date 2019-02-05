"""Module for IO operations"""

import numpy
import matplotlib.pyplot as plt
from matplotlib import cm

def plotContour(grid,var):
	"""
	Function to plot contours on a meshgrid

	Arguments
	---------

	grid : object of class Grid

	var  : integer
	     index of the plotting variable within grid data structure

	"""

	X,Y = numpy.meshgrid(grid.x_center,grid.y_center)
	X.transpose()
	Y.transpose()

	plt.figure()
	plt.contourf(X,Y,grid.data[:,:,var],cmap=cm.jet)
	plt.colorbar()
	plt.axis('equal')
	plt.show()
	
	return
