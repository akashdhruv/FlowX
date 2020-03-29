from qiskit.ignis.mitigation.measurement import complete_meas_cal, CompleteMeasFitter
from qiskit import IBMQ, Aer, BasicAer, execute
from qiskit.providers.ibmq import least_busy
from qiskit.tools.monitor import job_monitor
from qiskit.providers.aer.noise import NoiseModel

def calibrate_circuit_QASM(register,backend):
 
    provider = IBMQ.load_account()
    device_backend = provider.get_backend(backend)
    noise_model = NoiseModel.from_backend(device_backend)
    basis_gates = noise_model.basis_gates

    meas_calibs, state_lables = complete_meas_cal(qr=register, circlabel='mcal')

    device = Aer.get_backend('qasm_simulator')
    print("Running on device: ", device)

    job = execute(meas_calibs, backend=device, shots=1024, noise_model=noise_model, basis_gates=basis_gates)
    job_monitor(job, interval = 2)

    cal_results = job.result()

    meas_fitter = CompleteMeasFitter(cal_results, state_lables, circlabel='mcal')

    return meas_fitter, device, [noise_model, basis_gates]

def calibrate_circuit_IBMQ(register, backend):

    provider = IBMQ.load_account()

    meas_calibs, state_lables = complete_meas_cal(qr=register, circlabel='mcal')

    device = provider.get_backend(backend)
    print("Running on device: ", device)

    job = execute(meas_calibs, backend=device, shots=1024, max_credits=10)
    job_monitor(job, interval = 2)

    cal_results = job.result()

    meas_fitter = CompleteMeasFitter(cal_results, state_lables, circlabel='mcal')

    return meas_fitter, device, [None, None]
