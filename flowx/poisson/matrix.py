"""Routine to construct the matrix coefficient A of a linear system Ax = B"""


import numpy
from scipy import sparse
import time

start=time.time()
def construct_matrix(grid, user_bc):
    """Solve the Poisoon system using a Gauss-Seidel method
    
    Arguments
    ---------
    grid : Grid object
        Grid containing data.
    
    Returns 
    -------
    mtx : string
        Name of the coefficient matrix of the linear system.
    duration : float
        Time of the matrix construction 
    """
    
    #A-Diagonal
    A = numpy.ones(grid.nx*grid.ny)
    
    #B-Diagonal
    B = numpy.copy(A)
    B[grid.nx-1:((grid.nx*grid.ny)-1):grid.nx] = 0
    
    #D-Diagonal
    D = numpy.copy(A)
    D[grid.nx:((grid.nx*grid.ny)-1):grid.nx] = 0

    #C-Diagonal
    if(user_bc == 'dirichlet'):
        a = -6
        b = -5   
    else: 
        a = -2
        b = -3
    C1 = numpy.ones(grid.nx) * (b)
    C2 = numpy.copy(C1)
    C1[0:grid.nx:grid.nx-1] = a
    C2[1:grid.nx-1] = -4
    C3 = numpy.tile((C2), grid.nx-2)
    C = numpy.concatenate((C1, C3, C1), axis=0)

    #Matrix Construction
    data = numpy.array([A, B, C, D, A])
    diags = numpy.array([-grid.nx, -1, 0, 1, grid.nx])
    temp = sparse.spdiags(data, diags, grid.nx*grid.ny, grid.nx*grid.ny).toarray()
    temp = temp.astype(numpy.float64)

    #Matrix Compression
    mtx = sparse.csc_matrix(temp)
    
    end=time.time()
    duration=end-start
    return mtx, duration
