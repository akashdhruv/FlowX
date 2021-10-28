#!/usr/bin/env python
# coding: utf-8

# In[1]:


import flowx
from qiskit.visualization import plot_histogram
import numpy
flowx.__environment__


# In[2]:


nx, ny = 4, 4
qubits = numpy.log2(nx*ny)
xmin, xmax = -0.8, 0.5
ymin, ymax = -0.8, 0.5

center_vars   = ['lset']

particle_info = [dict(input='HDF5', file='sm_tshape.h5')]

quantum_info = dict(simulator = 'QASM', 
                    qubits = qubits, repeat = 1, circuit = 'grover', 
                    backend = 'ibmq_santiago', calibrate = False)

gridc, gridx, gridy, scalars, particles = flowx.domain_main(nx, ny, xmin, xmax, ymin, ymax, 
                                                            center_vars, particle_info=particle_info)

domain_data_struct = [gridc, gridx, gridy, scalars, particles]

grover = flowx.quantum_main(domain_data_struct, quantum_info=quantum_info)


# In[3]:


particles[0].rotate(45)
grover.setup_circuit()
grover.circuit.draw(output="mpl")


# In[4]:


flowx.io.plot_mesh_mapped(gridc, particles[0])


# In[5]:


results, answer = grover.run_circuit()
plot_histogram(answer['noisy'], bar_labels = False)


# In[6]:


plot_histogram([answer['noisy'], answer['mitigated']], legend=['noisy', 'mitigated'], bar_labels = False)

