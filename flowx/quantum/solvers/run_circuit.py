from qiskit import IBMQ, Aer, BasicAer, execute
from qiskit.providers.ibmq import least_busy
from qiskit.tools.monitor import job_monitor
from qiskit.providers.aer.noise import NoiseModel

def run_circuit_QASM(device, noise, circuit, quantum_register, classical_register):

    circuit.measure(quantum_register,classical_register)

    job = execute(circuit, backend=device, shots=1024, noise_model=noise[0], basis_gates=noise[1])
    job_monitor(job, interval = 2)

    results = job.result()
    answer = results.get_counts()

    return results, answer

def run_circuit_IBMQ(device, noise, circuit, quantum_register, classical_register):

    circuit.measure(quantum_register,classical_register)

    job = execute(circuit, backend=device, shots=1024, max_credits=10)
    job_monitor(job, interval = 2)

    results = job.result()
    answer = results.get_counts()

    return results, answer
