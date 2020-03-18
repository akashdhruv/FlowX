import numpy
import h5py

class Particles(object):
    """
    Class to store and advance particle data related to the immersed boundary 
    """

    def __init__(self,particle_info, xmin, xmax, ymin, ymax, scalars):
        """
        This constuctor allows user to initialize immersed boundary particles
        specific to their simulation.
        """

        self._scalars = scalars

        self.xmin = numpy.array([xmin, ymin]) 
        self.xmax = numpy.array([xmax, ymax])

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

        self.tracing_tol = 0.5*numpy.max(numpy.sqrt((self.x[2:,0]-self.x[1:-1,0])**2 + (self.x[2:,1]-self.x[1:-1,1])**2))

        # Procedure to find member variables
        #members = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("__")]

    def advance(self):
        """
        Subroutine to advance the particle data
        """
        self.x = self.xo + self._scalars.time*self.vel


    def offset(self, transform):
        """
        Subroutine to rotate particles
        """
        self.x = self.x + numpy.array(transform)

    def rotate(self, alpha):
        """
        Subroutine to rotate particles
        """
        rotation = numpy.zeros((2,2))

        rotation[0,:] = numpy.array([numpy.cos(alpha*numpy.pi/180),-numpy.sin(alpha*numpy.pi/180)])
        rotation[1,:] = numpy.array([numpy.sin(alpha*numpy.pi/180),numpy.cos(alpha*numpy.pi/180)])

        self.x = numpy.matmul(self.x,rotation)

    def reset(self):
        """
        Subroutine to reset particle positions to their original values
        """
        self.x = self.xo
