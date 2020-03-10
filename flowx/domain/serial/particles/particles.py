import numpy
import h5py

class Particles(object):
    """
    Class to store and advance particle data related to the immersed boundary 
    """

    def __init__(self,particle_info):
        """
        This constuctor allows user to initialize immersed boundary particles
        specific to their simulation.
        """

        particle_dict = h5py.File(particle_info['file'],'r')
        mesh = particle_dict['mesh']

        self.nnp = numpy.array(mesh['nnp'])[0]
 
        xp = numpy.array(mesh['x'])
        xp = numpy.reshape(xp, (self.nnp, 1))

        yp = numpy.array(mesh['y'])
        yp = numpy.reshape(yp, (self.nnp, 1))

        self.xo = numpy.concatenate((xp,yp), axis=1)      

        self.x = self.xo
        self.vel = numpy.zeros_like(self.xo)
        self.freq = numpy.zeros_like(self.xo)

        self.vel[:,0] = particle_info['vel'][0]
        self.vel[:,1] = particle_info['vel'][1]

        # Procedure to find member variables
        #members = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("__")]

    def advance(self,scalars):
        """
        Subroutine to advance the particle data
        """

        self.x = self.xo + scalars.time*self.vel
