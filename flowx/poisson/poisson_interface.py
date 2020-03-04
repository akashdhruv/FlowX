"""Interface for Poisson solver module"""

import abc

class poisson_interface(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, '__init__') and
                callable(subclass.__init__) and
                hasattr(subclass, 'solve_poisson') and 
                callable(subclass.solve_poisson) or 
                NotImplemented)


    @abc.abstractmethod
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
        raise NotImplementedError()

    @abc.abstractmethod
    def solve_poisson(self, grid):
        """ Subroutine to solve poisson equation

        Arguments
        ---------

        grid : object

             Grid object where the poisson equation needs to be solved

        """

        raise NotImplementedError()
