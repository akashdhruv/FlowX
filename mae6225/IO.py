import numpy
import matplotlib.pyplot as plt
from matplotlib import cm

def plotContour(grid,var):
	X,Y = numpy.meshgrid(grid.x_center,grid.y_center)
	X.transpose()
	Y.transpose()

	plt.figure()
	plt.contourf(X,Y,grid.data[:,:,var],cmap=cm.jet)
	plt.colorbar()
	plt.axis('equal')
	plt.show()
	
	return
