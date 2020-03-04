class Particles(object):
    """
    Class to store and advance particle data related to the immersed boundary 
    """

    def __init__(self,particle_info):
        """
        This constuctor allows user to initialize immersed boundary particles
        specific to their simulation.
        """
        
        self.xo, self.yo, self.ux, self.uy, self.omega = [0.0]*5        

    def advance(self, scalars):
        """
        Subroutine to advance the particle data
        """
