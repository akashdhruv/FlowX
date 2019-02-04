import numpy

def getAnalytical(grid,asol):

    for i in range(grid.nx):
	for j in range(grid.ny):
    		grid.data[asol,i,j] = numpy.cos(2*numpy.pi*grid.x_cell[i])*numpy.cos(2*numpy.pi*grid.y_cell[j])

def createRHS(grid,rvar):

def solveJacobi(grid,ivar,rvar):

def getError(grid,err,ivar,asol): 

def applyBC:	
