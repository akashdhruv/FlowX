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
    
    circuit.cu1(numpy.pi/4, register[0], register[3])
    circuit.cx(register[0], register[1])
    circuit.cu1(-numpy.pi/4, register[1], register[3])
    circuit.cx(register[0], register[1])
    circuit.cu1(numpy.pi/4, register[1], register[3])
    circuit.cx(register[1], register[2])
    circuit.cu1(-numpy.pi/4, register[2], register[3])
    circuit.cx(register[0], register[2])
    circuit.cu1(numpy.pi/4, register[2], register[3])
    circuit.cx(register[1], register[2])
    circuit.cu1(-numpy.pi/4, register[2], register[3])
    circuit.cx(register[0], register[2])
    circuit.cu1(numpy.pi/4, register[2], register[3])

def oracle_gate(circuit, register, particles):
    """
    Oracle gate 
    
    Arguments
    ---------
    
    circuit : Quantum Circuit
    register : Quantum register  
    particles : particles 
    """
   
    polygon = Polygon(numpy.ndarray.tolist(particles[0].x[1:]))
    grid_points = [Point(x) for x in numpy.ndarray.tolist(particles[1].x)]

    oracle_mask = []

    for grid_point, qubit in zip(grid_points, register): 
        if not grid_point.within(polygon): oracle_mask.append(qubit) 
   
    circuit.x(oracle_mask)
    controlled_Z_gate(circuit, register)
    circuit.x(oracle_mask)

def amplification_gate(circuit, register, particles):
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
