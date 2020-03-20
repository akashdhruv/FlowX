"""Implementation of the immersed boundary module"""

from flowx.quantum.quantum_interface import quantum_interface

class quantum_main(quantum_interface):

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

        from flowx.quantum.solvers.initialize import initialize_quantum_system
        from flowx.quantum.solvers.grover import oracle_gate, amplification_gate
        from flowx.quantum.solvers.run_circuit import run_circuit_QASM, run_circuit_IBMQ

        self._options = {'simulator' : 'QASM', 'qubits': 4, 'repeat' : 1, 'circuit' : 'grover'}
        self._simulators = {'QASM' : run_circuit_QASM, 'IBMQ' : run_circuit_IBMQ}

        self._gridc, self._gridx, self._gridy, self._scalars, self._particles = domain_data_struct
 
        if quantum_info:
            for key in quantum_info: self._options[key] = quantum_info[key]

        self.qubits = self._options['qubits']
        self.circuit, self.quantum_register, self.classical_register = initialize_quantum_system(self.qubits)

        if self._options['circuit'] is 'grover':
            self._gates = [oracle_gate, amplification_gate]*self._options['repeat']

        self._run_circuit = self._simulators[self._options['simulator']]

        return

    def setup_circuit(self):
        """
        """
       
        for gate in self._gates: gate(self.circuit, self.quantum_register, self._particles)

        return

    def run_circuit(self):
        """
        """

        results, answer = self._run_circuit(self.circuit, self.quantum_register, self.classical_register)

        return results, answer
