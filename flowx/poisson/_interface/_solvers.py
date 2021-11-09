"""Routine to solve the Poisson equation."""

import numpy
import scipy.sparse as sps
from scipy.sparse.linalg.dsolve import linsolve
from scipy.sparse import spdiags, csr_matrix
from scipy.sparse import linalg as sla

def solve_stub(grid, ivar, rvar, options):
    """Stub for poisson solver.

    Arguments
    ---------
    grid : Grid object
        Grid containing data.
    ivar : string
        Name of the grid variable of the numerical solution.
    rvar : string
        Name of the grid variable of the right-hand side.
    options : dictionary

    Returns
    -------
    ites: None
        Number of iterations computed.
    residual: None
        Final residual.
    """

    return None, None

def solve_cg(grid, ivar, rvar, options):
    """Solve the Poisson system using a conjugate-gradient method.

    Arguments
    ---------
    grid : Grid object
        Grid containing data.
    ivar : string
        Name of the grid variable of the numerical solution.
    rvar : string
        Name of the grid variable of the right-hand side.

    options: dictionary

    Returns
    -------
    ites: integer
        Number of iterations computed.
    residual: float
        Final residual.
    """

    verbose = options['verbose']
    maxiter = options['maxiter']
    tol = options['tol']

    def A(p):
        return ((p[1:-1, :-2] - 2 * p[1:-1, 1:-1] + p[1:-1, 2:]) / dx**2 +
                (p[:-2, 1:-1] - 2 * p[1:-1, 1:-1] + p[2:, 1:-1]) / dy**2)

    def fill_guard_cells_neumann(x, bc_val, dx, dy):
        x[:, 0] = bc_val * dx + x[:, 1]
        x[:, -1] = bc_val * dx + x[:, -2]
        x[0, :] = bc_val * dy + x[1, :]
        x[-1, :] = bc_val * dy + x[-2, :]

    p = grid[ivar][0,0,:,:]  # initial guess
    b = grid[rvar][0,0,:,:]  # RHS of the system
    dx, dy = grid.dx, grid.dy  # cell widths

    r = b[1:-1, 1:-1] - A(p)  # initial residuals
    rk_norm = numpy.sum(r * r)  # inner product
    d = numpy.zeros_like(p)  # search direction
    d[1:-1, 1:-1] = r  # set direction to initial residual

    # Check for Neumann boundary conditions.
    bc_type, bc_val = grid.bc_type[ivar][0], grid.bc_val[ivar][0]
    if bc_type == 'neumann':
        # Apply Neumann boundary conditions to search direction.
        fill_guard_cells_neumann(d, bc_val, dx, dy)

    ites = 0  # iteration index
    res = rk_norm  # initial residual
    while ites < maxiter and res > tol:
        Ad = A(d)
        alpha = rk_norm / numpy.sum(d[1:-1, 1:-1] * Ad)  # step size
        p[1:-1, 1:-1] += alpha * d[1:-1, 1:-1]  # update solution
        r -= alpha * Ad  # update residuals
        r_norm = numpy.sum(r * r)  # inner product
        beta = r_norm / rk_norm
        rk_norm = r_norm
        d[1:-1, 1:-1] = r + beta * d[1:-1, 1:-1]  # update search direction
        if bc_type == 'neumann':
            fill_guard_cells_neumann(d, bc_val, dx, dy)
        res = r_norm
        ites += 1

    grid.fill_guard_cells(ivar)

    if verbose:
        print('CG method:')
        if ites == maxiter:
            print('Warning: maximum number of iterations reached!')
        print('- Number of iterations: {}'.format(ites))
        print('- Final residual: {}'.format(res))

    return ites, res

def solve_direct(grid, ivar, rvar, options):

    """Solve the Poisson system using a direct solver from the scipy library.

    Arguments
    ---------
    grid : Grid object
        Grid containing data.
    ivar : string
        Name of the grid variable of the numerical solution.
    rvar : string
        Name of the grid variable of the right-hand side.

    options : dictionary

    """
    verbose = options['verbose']
    matrix = options['matrix']

    dx, dy = grid.dx, grid.dy
    nx, ny = grid.nx, grid.ny

    rhs = grid[rvar][0,0,:,:]
    phi = grid[ivar][0,0,:,:]

    sol = linsolve.spsolve(matrix, rhs[1:-1, 1:-1].flatten())
    residual = numpy.linalg.norm(matrix * sol - rhs[1:-1, 1:-1].flatten())

    phi[1:-1,1:-1] = numpy.reshape(sol,(ny,nx))
    grid.fill_guard_cells(ivar)

    if verbose:
        print('Direct Solver:')
        print('- Final residual: {}'.format(residual))

    return None, residual

def solve_jacobi(grid, ivar, rvar, options):
    """Solve the Poisson system using a Jacobi method.

    Arguments
    ---------
    grid : Grid object
        Grid containing data.
    ivar : string
        Name of the grid variable of the numerical solution.
    rvar : string
        Name of the grid variable of the right-hand side.
    options : dictionary

    Returns
    -------
    ites: integer
        Number of iterations computed.
    residual: float
        Final residual.
    """

    maxiter = options['maxiter']
    tol = options['tol']
    verbose = options['verbose']

    phi = grid[ivar][0,0,:,:]
    b = grid[rvar][0,0,:,:]
    dx, dy = grid.dx, grid.dy

    ites = 0
    residual = tol + 1.0
    while ites < maxiter and residual > tol:
        phi_old = numpy.copy(phi)  # previous solution
        phi[1:-1, 1:-1] = (((phi_old[:-2, 1:-1] +
                             phi_old[2:,  1:-1]) * dy**2 +
                            (phi_old[1:-1, :-2] +
                             phi_old[1:-1,  2:]) * dx**2 -
                            b[1:-1, 1:-1] * dx**2 * dy**2) /
                           (2 * (dx**2 + dy**2)))

        grid.fill_guard_cells(ivar)

        residual = (numpy.sqrt(numpy.sum((phi - phi_old)**2) /
                    ((grid.nx + 2) * (grid.ny + 2))))
        ites += 1

    if verbose:
        print('Jacobi method:')
        if ites == maxiter:
            print('Warning: maximum number of iterations reached!')
        print('- Number of iterations: {}'.format(ites))
        print('- Final residual: {}'.format(residual))

    return ites, residual

def solve_superlu(grid, ivar, rvar, options):
    """Solve the Poisson system using a direct solver from the scipy library.

    Arguments
    ---------
    grid : Grid object
        Grid containing data.
    ivar : string
        Name of the grid variable of the numerical solution.
    rvar : string
        Name of the grid variable of the right-hand side.

    options : dictionary

    """

    verbose = options['verbose']
    matrix = options['matrix']
    lu = options['lu']

    dx, dy = grid.dx, grid.dy
    nx, ny = grid.nx, grid.ny

    rhs = grid[rvar][0,0,:,:]
    phi = grid[ivar][0,0,:,:]

    sol = lu.solve(rhs[1:-1, 1:-1].flatten())
    residual = numpy.linalg.norm(matrix * sol - rhs[1:-1, 1:-1].flatten())

    phi[1:-1,1:-1] = numpy.reshape(sol,(ny,nx))
    grid.fill_guard_cells(ivar)

    if verbose:
        print('LU Decomposition:')
        print('- Final residual: {}'.format(residual))

    return None, residual
