"""Module with iterative solvers for the Poisson system."""

import numpy


def solve_jacobi(grid, ivar, rvar, maxiter=3000, tol=1e-9):
    """Solve the Poisson system using a Jacobi method.

    Arguments
    ---------
    grid : Grid object
        Grid containing data.
    ivar : string
        Name of the grid variable of the numerical solution.
    rvar : string
        Name of the grid variable of the right-hand side.
    maxiter : integer, optional
        Maximum number of iterations;
        default: 3000
    tol : float, optional
        Exit-criterion tolerance;
        default: 1e-9

    Returns
    -------
    ites: integer
        Number of iterations computed.
    residual: float
        Final residual.

    """
    i_ivar, i_rvar = grid.get_variable_indices([ivar, rvar])
    dx, dy = grid.dx, grid.dy

    ites = 0
    residual = tol + 1.0
    while ites < maxiter and residual > tol:
        phi_old = numpy.copy(grid.data[:, :, i_ivar])  # previous solution
        grid.data[1:-1, 1:-1, i_ivar] = (((phi_old[1:-1, :-2] +
                                           phi_old[1:-1, 2:]) * dy**2 +
                                          (phi_old[:-2, 1:-1] +
                                           phi_old[2:, 1:-1]) * dx**2 -
                                          grid.data[1:-1, 1:-1, i_rvar] *
                                          dx**2 * dy**2) /
                                         (2.0 * (dx**2 + dy**2)))

        grid.fill_guard_cells(ivar)

        residual = (numpy.sqrt(numpy.sum((grid.data[:, :, i_ivar] -
                                          phi_old)**2) /
                    ((grid.nx + 2) * (grid.ny + 2))))
        ites += 1

    if ites == maxiter:
        print('Warning: maximum number of iterations reached!')

    return ites, residual
