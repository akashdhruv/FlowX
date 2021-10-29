"""Implementation of the quantum module"""

from . import _interface

class Quantum(object):

    def __init__(self, domain_data_struct=[None]*5, quantum_vars=[None]*2, quantum_info=None):

        """
        Constructor for the quantum unit

        Arguments
        ---------

        domain_data_struct : object list
              [gridc, gridx, gridy, scalars, particles]

        quantum_vars : list
                List of string for field variables required by quantum unit
               
        qauntum_info : Dictionary of keyword arguments

        """
        #----------------Create images of other units and variables--------------------------------
        self._gridc, self._gridx, self._gridy, self._scalars, self._particles = domain_data_struct
        self._ibmf, self._velc = quantum_vars
 
        #----------------Setup default parameters-------------------------------------------------
        self._options = {'simulator' : 'QASM', \
                         'qubits': 4, \
                         'repeat' : 1, \
                         'circuit' : 'grover', \
                         'backend': 'ibmq_london', \
                         'calibrate' : False}

        self._simulators = {'QASM' : _interface.run_circuit_QASM, \
                            'IBMQ' : _interface.run_circuit_IBMQ}

        self._calibrators = {'QASM' : _interface.calibrate_circuit_QASM, \
                             'IBMQ' : _interface.calibrate_circuit_IBMQ}

        #----------------Read user parameters--------------------------------------------
        if quantum_info:
            for key in quantum_info: self._options[key] = quantum_info[key]

        #----------------Setup current unit-----------------------------------------------

        if None in domain_data_struct: #or None in quantum_vars:
            raise ValueError('Quantum unit cannot be setup. One or more parameter is missing')

        else:
            self.qubits = self._options['qubits']

            self.circuit, self.quantum_register, self.classical_register = \
            _interface.initialize_quantum_system(self.qubits)

            self._calibrate_circuit = self._calibrators[self._options['simulator']]
            self._run_circuit = self._simulators[self._options['simulator']]

            self.fitter, self.device, self.noise = \
            self._calibrate_circuit(self.quantum_register, self._options['backend'], 
                                                           self._options['calibrate'])

            if self._options['circuit'] == 'grover':
                self._gates = [_interface.grover.oracle_gate, 
                               _interface.grover.amplification_gate]*self._options['repeat']

        return

    def setup_circuit(self):
        """
        """
       
        for gate in self._gates: gate(self.circuit, self.quantum_register, self._particles, self._gridc)

        return

    def run_circuit(self):
        """
        """

        results, answer = self._run_circuit(self.device, self.noise, self.fitter, self.circuit, self.quantum_register, self.classical_register)

        return results, answer
