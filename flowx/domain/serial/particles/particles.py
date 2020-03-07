import numpy

class Particles(object):
    """
    Class to store and advance particle data related to the immersed boundary 
    """

    def __init__(self,particle_info):
        """
        This constuctor allows user to initialize immersed boundary particles
        specific to their simulation.
        """

        self.xo, self.radius, self.vel, self.freq = [0.0]*4

        for key, value in particle_info.items(): setattr(self, key, numpy.array(value))

        self.x = self.xo

        #members = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("__")]

        #for member in members:
        #    setattr(self, member, numpy.array(getattr(self, member))) 

    def advance(self,scalars):
        """
        Subroutine to advance the particle data
        """

        if self.freq:
            self.x = self.xo + numpy.sin(2.*numpy.pi*self.freq*scalars.time) 
            self.vel = 2.*numpy.pi*self.freq*numpy.cos(2.*numpy.pi*self.freq*scalars.time)

        else:
            self.x = self.xo + scalars.time*self.vel
