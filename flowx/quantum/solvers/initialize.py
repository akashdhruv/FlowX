from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister

def initialize_quantum_system(qubits):
    """
    Initialize the quantum system
    
    Arguments
    ---------
    
    qubits : Number of quantum bits
    """
    
    quantum_register = QuantumRegister(qubits,'q')
    classical_register = ClassicalRegister(qubits,'c')
    calibration_register = QuantumRegister(qubits, 'q')
    
    circuit = QuantumCircuit(quantum_register, classical_register)   
    circuit.h(quantum_register)
  
    return circuit, quantum_register, classical_register, calibration_register
