import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg.dsolve import linsolve


def solve_direct(grid, ivar, rvar, verbose=False):

    """Solve the Poisson system using a direct solver from the scipy library.

    Arguments
    ---------
    grid : Grid object
        Grid containing data.
    ivar : string
        Name of the grid variable of the numerical solution.
    rvar : string
        Name of the grid variable of the right-hand side.

    Returns
    -------
    verbose : bool, optional
        Set True to display residual information;
        default: False.

    """

    phi = grid.get_values(ivar)
    b = grid.get_values(rvar)
    dx, dy = grid.dx, grid.dy
    nx, ny = grid.nx, grid.ny

    matrix_length = nx*ny

    mtx = sps.lil_matrix((matrix_length, matrix_length), dtype=np.float64)

    counter = 0

    for i in range(1,nx+1):
        for j in range(1,ny+1):

            coeff = -8.0
    
            if(j > 1):
                mtx[counter,counter-1] = 1.0
                coeff = coeff + 1.0

            if(j < ny):
                mtx[counter,counter+1] = 1.0
                coeff = coeff + 1.0

            if(i > 1):
                mtx[counter,counter-ny] = 1.0
                coeff = coeff + 1.0

            if(i < nx):
                mtx[counter,counter+ny] = 1.0
                coeff = coeff + 1.0

            mtx[counter,counter] = coeff

            counter = counter + 1
        

    mtx = mtx/(dx**2)
    mtx = mtx.tocsr()
    rhs = b[1:-1, 1:-1].flatten()
    sol = linsolve.spsolve(mtx, rhs)

    residual = np.linalg.norm(mtx * sol - rhs)

    phi[1:-1,1:-1] = np.reshape(sol,(nx,ny))
    grid.fill_guard_cells(ivar)

    if verbose:
        print('Direct Solver:')
        print('- Final residual: {}'.format(residual))

    return residual
