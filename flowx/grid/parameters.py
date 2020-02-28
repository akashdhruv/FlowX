"""Module for simulation constants"""

import numpy

class Parameters(object):

    def __init__(self, tmax = 0.0, dt = 1.0, Re = 1.0):
        
        self.tmax = tmax
        self.dt = dt
        self.Re = Re
        self.to = 0.0
        self.nstep = 0
        self.time = 0

    def advance(self):
        self.time += self.dt
        self.nstep += 1
