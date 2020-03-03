"""Interface for incompressible Navier Stokes module"""

import abc

class ins_interface(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, '__init__') and
                callable(subclass.__init__) and
                hasattr(subclass, 'advance') and 
                callable(subclass.advance) or 
                NotImplemented)


    @abc.abstractmethod
    def __init__(self, ins_vars, **kwargs):

        """
        Constructor for the ins unit

        Arguments
        ---------

        ins_vars : list
                List of string for field variables required by ins unit
               
                ins_vars[0] --> velocity
                ins_vars[1] --> RHS from the previous time step
                ins_vars[2] --> divergence
                ins_vars[3] --> pressure


        **kwargs : Dictionary of keyword arguments

        'time_stepping' keyword refers to the time advancement scheme to be used

        kwargs['time_stepping'] = 'ab2' --> default
                                = 'euler'

        """
        raise NotImplementedError()

    @abc.abstractmethod
    def advance(self, poisson, imbound, gridc, gridx, gridy, scalars, particles):

        """
        Subroutine for the fractional step explicit time advancement of Navier Stokes equations
 
        Arguments
        ---------
        poisson : object
            Object for the poisson solver

        imbound : object
            Object for the immersed boundary unit

        gridc : object
          Grid object for cell centered variables

        gridx : object
          Grid object for x-face variables

        gridy : object
          Grid object for y-face variables

        scalars: object
           Scalars object to access time-step and Reynold number

        particles: object
           Object containing immersed boundary information
        """

        raise NotImplementedError()
