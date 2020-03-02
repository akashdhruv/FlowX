"""Routine to solve the Poisson equation."""

from flowx.poisson.solvers.jacobi import solve_jacobi

def poisson(grid, ivar, rvar, options=None):

    ites,residual = solve_jacobi(grid, ivar, rvar, maxiter=3000, tol=1e-9, verbose=False)

    return ites,residual
