class Particles(object):
    """
    Class to store and advance particle data related to the immersed boundary 
    """

    def __init__(self,particle_info):
        """
        This constuctor allows user to initialize immersed boundary particles
        specific to their simulation.
        """

        self.xo, self.yo, self.radius, self.velx, self.vely = [0.0]*5

        if 'xo' in particle_info: self.xo = particle_info['xo']
        if 'yo' in particle_info: self.yo = particle_info['yo']
        if 'radius' in particle_info: self.radius = particle_info['radius']
        if 'velx' in particle_info: self.velx = particle_info['velx']
        if 'vely' in particle_info: self.vely = particle_info['vely']


        self.x = self.xo
        self.y = self.yo
      
    def advance(self, scalars):
        """
        Subroutine to advance the particle data
        """

        self.x = self.xo + scalars.variable['time']*self.velx
        self.y = self.yo + scalars.variable['time']*self.vely
        
