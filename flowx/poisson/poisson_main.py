"""Routine to solve the Poisson equation."""

from flowx.poisson.solvers.serial.jacobi import solve_serial_jacobi
from flowx.poisson.solvers.serial.cg import solve_serial_cg
from flowx.poisson.solvers.parallel.jacobi import solve_parallel_jacobi
from flowx.poisson.solvers.parallel.cg import solve_parallel_cg

def solve_poisson(grid, ivar, rvar, **kwargs):

    _solver_type = 'serial_cg'

    if 'poisson_solver' in kwargs: _solver_type = kwargs.get('poisson_solver')

    if _solver_type is 'serial_jacobi':
        solve_poisson = solve_serial_jacobi

    elif _solver_type is 'serial_cg':
        solve_poisson = solve_serial_cg 

    elif _solver_type is 'parallel_jacobi':
        solve_poisson = solve_parallel_jacobi 

    elif _solver_type is 'parallel_cg':
        solve_poisson = solve_parallel_cg

    ites,residual = solve_poisson(grid, ivar, rvar, **kwargs)

    return ites,residual
