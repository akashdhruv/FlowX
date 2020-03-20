"""Interface for quantum computing module"""

import abc

class quantum_interface(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, '__init__') and
                callable(subclass.__init__) and
                hasattr(subclass, 'setup_circuit') and
                callable(subclass.setup_circuit) and
                hasattr(subclass, 'run_circuit') and
                callable(subclass.run_circuit) or
                NotImplemented)

    @abc.abstractmethod
    def __init__(self, domain_data_struct=[None]*5, quantum_vars=None, quantum_info=None):

        """
        Constructor for the imbound unit

        Arguments
        ---------
        domain_data_struct : object list
               [gridc, gridx, gridy, scalars, particles]

        quantum_vars : list
                List of string for field variables required by ins unit

        quantum_info : Dictionary of keyword arguments

        """
        raise NotImplementedError()

    @abc.abstractmethod
    def setup_circuit(self):
        """
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def run_circuit(self):
        """
        """
        raise NotImplementedError()
