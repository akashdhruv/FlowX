"""Interface for Poisson solver module"""

from flowx.poisson.poisson_interface import poisson_interface

class poisson_main(poisson_interface):

    def __init__(self, poisson_vars=None, poisson_info=None):

        """
        Constructor for the Poisson unit

        Arguments
        ---------

        poisson_vars : list
                List of string for field variables required by poisson unit
               
                poisson_vars[0] --> Phi (numerical solution)
                poisson_vars[1] --> RHS

        poisson_info : Dictionary of keyword arguments

        'poisson_solver' keyword refers to the type of solver to be used
        poisson_info['poisson_solver'] = 'serial_cg' --> default
                                       = 'serial_jacobi'

        poisson_info['maxiter'] = maximum number of iterations --> default 2000
       
        poisson_info['tol']  = minimum tolerance of the residuals --> default 1e-9
 
        poisson_info['verbose'] = bool to displacy poisson stats or not --> default False

        """

        from flowx.poisson.solvers.serial.jacobi import solve_serial_jacobi
        from flowx.poisson.solvers.serial.cg import solve_serial_cg

        self._ivar = 'stub'
        self._rvar = 'stub'

        self._solver_type = 'serial_cg'
        self._maxiter = 2000
        self._tol = 1e-9
        self._verbose = False

        if poisson_info:
            if 'poisson_solver' in poisson_info: self._solver_type = poisson_info['poisson_solver']
            if 'maxiter' in poisson_info: self._maxiter = poisson_info['maxiter']
            if 'tol' in poisson_info: self._tol = poisson_info['tol']
            if 'verbose' in poisson_info: self._verbose = poisson_info['verbose']

        if self._solver_type is 'serial_cg':
            self._solve_poisson = solve_serial_cg
        elif self._solver_type is 'serial_jacobi':
            self._solve_poisson = solve_serial_jacobi

        if poisson_vars:
            self._ivar = poisson_vars[0]
            self._rvar = poisson_vars[1]

        else:
            print('Warning: Poisson unit is a stub, any call to its methods will result in an error.') 

        return

    def solve_poisson(self, grid):
        """ Subroutine to solve poisson equation

        Arguments
        ---------

        grid : object

             Grid object where the poisson equation needs to be solved

        """

        ites, residual = self._solve_poisson(grid, self._ivar, self._rvar, self._maxiter, self._tol, self._verbose)

        return ites, residual
