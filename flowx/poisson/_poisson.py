"""Interface for Poisson solver module"""

from . import _interface

class Poisson(object):

    def __init__(self, grid=None, poisson_vars=[None]*2, poisson_info=None):

        """
        Constructor for the Poisson unit

        Arguments
        ---------

        grid : object

             Grid object where the poisson equation needs to be solved

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
        #---------------------Create images of other units and objects----------------
        self._grid = grid
        self._ivar, self._rvar = poisson_vars

        #--------------------Set default parameters---------------------------------
        self._options = {'poisson_solver' : 'superlu', \
                         'maxiter': 2000, \
                         'tol' : 1e-9, \
                         'verbose' : False}

        self._serial_iterative_solvers = {'cg' : _interface.solve_cg, \
                                          'jacobi': _interface.solve_jacobi}

        self._serial_direct_solvers = {'direct' : _interface.solve_direct, \
                                       'superlu' : _interface.solve_superlu}

        #----------------------Read user parameters------------------------------------
        if poisson_info: 
            for key in poisson_info: self._options[key] = poisson_info[key]

        #----------------------Setup current unit with default/user parameters-----------

        if not grid or None in poisson_vars:
            self._solve = _interface.solveStub
            print('Warning: Poisson unit is a stub') 

        else:
            self._solve = {**self._serial_iterative_solvers, 
                           **self._serial_direct_solvers}[self._options['poisson_solver']]

            self._options['lu'], self._options['mtx'] = _interface.build_sparse_matrix(self._grid, self._ivar)

        return

    def solve(self):
        """ Subroutine to solve poisson equation

        """

        ites, residual = self._solve(self._grid, self._ivar, self._rvar, self._options)

        return ites, residual
