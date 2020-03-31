import numpy
from shapely.geometry import Point, Polygon

def controlled_Z_gate(circuit, register):
    """
    Apply a triple controlled Z gate to qubit 3 
    
    Arguments
    ---------
    
    circuit : Quantum Circuit
    register : Quantum register
    """
    n_qubits = len(register)

    circuit.cu1(numpy.pi/4, register[0], register[n_qubits-1])

    for qubit in range(n_qubits-1):
        
        if qubit < n_qubits-2:
            circuit.cx(register[qubit], register[qubit+1])
            circuit.cu1(-numpy.pi/4, register[qubit+1], register[n_qubits-1])

            circuit.cx(register[0], register[qubit+1])
            circuit.cu1(numpy.pi/4, register[qubit+1], register[n_qubits-1])

        else:
            circuit.cx(register[qubit-1], register[qubit])
            circuit.cu1(-numpy.pi/4, register[qubit], register[n_qubits-1])

            circuit.cx(register[0], register[qubit])
            circuit.cu1(numpy.pi/4, register[qubit], register[n_qubits-1])

def oracle_gate(circuit, register, particles, grid):
    """
    Oracle gate 
    
    Arguments
    ---------
    
    circuit : Quantum Circuit
    register : Quantum register  
    particles : particles 
    """
   
    polygon = Polygon(numpy.ndarray.tolist(particles[0].x[1:]))

    nx, ny = grid.nx, grid.ny
    x, y = grid.x, grid.y

    for i in range(1,nx+1):
        for j in range(1,ny+1):   

            grid_point = Point(x[i], y[j])

            x_bin = "{0:b}".format(i-1)
            y_bin = "{0:b}".format(j-1)

            x_append = '0'*(int(len(register)/2) - len(x_bin))
            y_append = '0'*(int(len(register)/2) - len(y_bin))

            x_bin = "".join([x_append, x_bin])
            y_bin = "".join([y_append, y_bin])

            oracle_list = list("".join([y_bin, x_bin]))
            list.reverse(oracle_list)

            oracle_mask = []

            if(grid_point.within(polygon)):

                for oracle_char, qubit in zip(oracle_list, register):
                    if(oracle_char is '0'): oracle_mask.append(qubit)

                if(oracle_mask): circuit.x(oracle_mask)
                controlled_Z_gate(circuit, register)
                if(oracle_mask): circuit.x(oracle_mask)

def amplification_gate(circuit, register, particles, grid):
    """
    Amplification_gate
    
    Arguments
    ---------
    
    circuit : Quantum Circuit
    register : Quantum register
    particles : particles        
    """
    
    circuit.h(register)
    circuit.x(register)
    
    controlled_Z_gate(circuit, register)
    
    circuit.x(register)
    circuit.h(register)
