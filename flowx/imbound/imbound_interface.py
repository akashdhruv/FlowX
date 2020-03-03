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
                callable(subclass.map_to_grid) or
                NotImplemented)


    @abc.abstractmethod
    def __init__(self, imbound_vars=None, **kwargs):

        """
        Constructor for the imbound unit

        Arguments
        ---------

        imbound_vars : list
                List of string for field variables required by ins unit
               
                imbound_vars[0] --> indicator variable for immersed boundary
                imbound_vars[1] --> velocity variable

        **kwargs : Dictionary of keyword arguments

        'with_ib' - keyword to indicate if immersed boundary is present or not
                  - True if present
                  - False if not present

        kwargs['with_ib'] = False --> default
                          = True

        Error generated if kwargs['with_ib'] is True and imbound_vars is None

        """
        raise NotImplementedError()

    @abc.abstractmethod
    def map_to_grid(self, gridx, gridy, particles):
        """
        Subroutine to map immersed boundary on grid
 
        Arguments
        ---------
        gridx : object
          Grid object for x-face variables

        gridy : object
          Grid object for y-face variables

        particles: object
           Object containing immersed boundary information
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def force_flow(self, gridx, gridy, scalars, particles):

        """
        Subroutine to compute immersed boundary forces
 
        Arguments
        ---------
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
