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
    def __init__(self, poisson=None, imbound=None, domain_data_struct=[None]*5, ins_vars=[None]*5, ins_info=None):

        """
        Constructor for the ins unit

        Arguments
        ---------

        poisson : object
            Object for the poisson solver

        imbound : object
            Object for the immersed boundary unit

        domain_data_struct : object list
          [gridc, gridx, gridy, scalars, particles]


        ins_vars : list
                List of string for field variables required by ins unit
               
                ins_vars[0] --> velocity
                ins_vars[1] --> RHS from the previous time step
                ins_vars[2] --> divergence
                ins_vars[3] --> pressure
                ins_vars[4] --> delta pressure

        ins_info : Dictionary of keyword arguments

        'time_stepping' keyword refers to the time advancement scheme to be used

        ins_info['time_stepping'] = 'ab2' --> default
                                  = 'euler'

        """
        raise NotImplementedError()

    @abc.abstractmethod
    def advance(self):

        """
        Subroutine for the fractional step explicit time advancement of Navier Stokes equations
 
        """

        raise NotImplementedError()
