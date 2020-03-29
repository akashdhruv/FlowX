from qiskit import IBMQ, Aer, BasicAer, execute
from qiskit.providers.ibmq import least_busy
from qiskit.tools.monitor import job_monitor
from qiskit.providers.aer.noise import NoiseModel

def run_circuit_statevector(circuit, quantum_register, classical_register):
    backend_sim = Aer.get_backend('statevector_simulator')
    job_sim = execute(grover_circuit, backend_sim)
    statevec = job_sim.result().get_statevector()
    print(statevec)

    return None, None

def run_circuit_QASM(circuit, quantum_register, classical_register):

    circuit.measure(quantum_register,classical_register)

    backend = BasicAer.get_backend('qasm_simulator')
    shots = 1024
    results = execute(circuit, backend=backend, shots=shots).result()
    answer = results.get_counts()

    return results, answer

def run_circuit_IBMQ(circuit, quantum_register, classical_register):

    circuit.measure(quantum_register,classical_register)

    provider = IBMQ.load_account()
    device = least_busy(provider.backends(simulator=False))
    print("Running on current least busy device: ", device)

    job = execute(circuit, backend=device, shots=1024, max_credits=10)
    job_monitor(job, interval = 2)

    results = job.result()
    answer = results.get_counts(circuit)

    return results, answer
