"""Module with iterative solvers for the Poisson system."""

import numpy
 

def solve_SOR(grid, ivar, rvar, maxiter=3000, tol=1e-9, verbose=False):
    """Solve the Poisson system using a Successive Over-Relaxation method.

    Arguments
    ---------
    grid : Grid object
        Grid containing data.
    ivar : string
        Name of the grid variable of the numerical solution.
    rvar : string
        0 Name of the grid variable of the right-hand side.
    maxiter : integer, optional
        Maximum number of iterations;
        default: 3000
    tol : float, optional
        Exit-criterion tolerance;
        default: 1e-9
    omega : float, optional
        relaxation parameter.

    Returns
    -------
    ites: integer
        Number of iterations computed.
    residual: float
        Final residual.
    verbose : bool, optional
        Set True to display convergence information;
        default: False.

    """
    i_ivar, i_rvar = grid.get_variable_indices(ivar, rvar)
    nx, ny = grid.nx, grid.ny
    dx, dy = grid.dx, grid.dy

    ites = 0
    residual = tol + 1.0
    omega = 1.50
    
    while ites < maxiter and residual > tol:
        phi_old = numpy.copy(grid.data[:, :, i_ivar])  # previous solution
       
        for j in range(1,ny+1):
            for i in range(1,nx+1):
                grid.data[i, j, i_ivar] = omega*(((grid.data[i, j-1, i_ivar] +
                                                   grid.data[i, j+1, i_ivar]) * dy**2 +
                                                  (grid.data[i-1, j, i_ivar] +
                                                   grid.data[i+1, j, i_ivar]) * dx**2 -
                                                  grid.data[i, j, i_rvar] *
                                                  dx**2 * dy**2) /
                                                 (2.0 * (dx**2 + dy**2))) +(1.0 - omega) * grid.data[i, j, i_ivar]

        grid.fill_guard_cells(ivar)

        residual = (numpy.sqrt(numpy.sum((grid.data[:, :, i_ivar] -
                                          phi_old)**2) /
                    ((grid.nx + 2) * (grid.ny + 2))))
        ites += 1

    if verbose:
        print('Successive Over-Relaxation method:')
        if ites == maxiter:
            print('Warning: maximum number of iterations reached!')
        print('- Number of iterations: {}'.format(ites))
        print('- Final residual: {}'.format(residual))

    return ites, residual

