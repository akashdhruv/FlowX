"""Interface for immersed boundary module"""

import abc

class imbound_interface(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, '__init__') and
                callable(subclass.__init__) and
                hasattr(subclass, 'force_flow') and 
                callable(subclass.force_flow) and
                hasattr(subclass, 'map_to_grid') and
                callable(subclass.map_to_grid) and
                hasattr(subclass, 'advect') and
                callable(subclass.advect) or
                NotImplemented)


    @abc.abstractmethod
    def __init__(self, domain_data_struct=[None]*5, imbound_vars=[None]*4, imbound_info=None):

        """
        Constructor for the imbound unit

        Arguments
        ---------
        domain_data_struct : object list
               [gridc, gridx, gridy, scalars, particles]

        imbound_vars : list
                List of string for field variables required by ins unit
               
                imbound_vars[0] --> indicator variable for immersed boundary
                imbound_vars[1] --> velocity variable
                imbound_vars[2] --> dynamic x grid variable
                imbound_vars[3] --> dynamic y grid variable

        imbound_info : Dictionary of keyword arguments

        'with_ib' - keyword to indicate if immersed boundary is present or not
                  - True if present
                  - False if not present

        imbound_info['with_ib'] = False --> default
                                = True

        Error generated if kwargs['with_ib'] is True and imbound_vars is empty

        """
        raise NotImplementedError()

    @abc.abstractmethod
    def map_to_grid(self):
        """
        Subroutine to map immersed boundary on grid
 
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def force_flow(self):

        """
        Subroutine to compute immersed boundary forces
 
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def advect(self):

        """
        Subroutine to advect immersed boundary
 
        """
        raise NotImplementedError()
