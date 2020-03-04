"""Interface for Poisson solver module"""

from flowx.poisson.poisson_interface import poisson_interface

class poisson_main(poisson_interface):

    def __init__(self, poisson_vars, **kwargs):

        """
        Constructor for the Poisson unit

        Arguments
        ---------

        poisson_vars : list
                List of string for field variables required by poisson unit
               
                poisson_vars[0] --> Phi (numerical solution)
                poisson_vars[1] --> RHS

        **kwargs : Dictionary of keyword arguments

        'solver_type' keyword refers to the type of solver to be used
        kwargs['solver_type'] = 'serial_cg' --> default
                              = 'serial_jacobi'

        kwargs['maxiter'] = maximum number of iterations --> default 2000
       
        kwargs['tol']  = minimum tolerance of the residuals --> default 1e-9
 
        kwargs['verbose'] = bool to displacy poisson stats or not --> default False

        """

        from flowx.poisson.solvers.serial.jacobi import solve_serial_jacobi
        from flowx.poisson.solvers.serial.cg import solve_serial_cg

        self._ivar = poisson_vars[0]
        self._rvar = poisson_vars[1]

        self._solver_type = 'serial_cg'
        self._maxiter = 2000
        self._tol = 1e-9
        self._verbose = False

        if 'solver_type' in kwargs: self._solver_type = kwargs['solver_type']
        if 'maxiter' in kwargs: self._maxiter = kwargs['maxiter']
        if 'tol' in kwargs: self._tol = kwargs['tol']
        if 'verbose' in kwargs: self._verbose = kwargs['verbose']          


        if self._solver_type is 'serial_cg':
            self._solve_poisson = solve_serial_cg
        elif self._solver_type is 'serial_jacobi':
            self._solve_poisson = solve_serial_jacobi

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
