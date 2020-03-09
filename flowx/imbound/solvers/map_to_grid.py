import numpy

def map_to_grid_stub(gridx, gridy, particles, ibmf):

    """
    Stub subroutine to compute IB mapping on grid
 
    Arguments
    ---------
    gridx : object
      Grid object for x-face variables

    gridy : object
      Grid object for y-face variables

    particles: object
       Object containing immersed boundary information

    ibmf : string for forcing variable

    """

    return

def map_to_grid_levelset(gridx, gridy, particles, ibmf):

    """
    Subroutine to compute IB mapping on grid using the level set function
 
    Arguments
    ---------
    gridx : object
      Grid object for x-face variables

    gridy : object
      Grid object for y-face variables

    particles: object
       Object containing immersed boundary information

    ibmf : string for forcing variable

    """

    Xfx, Yfx = numpy.meshgrid(gridx.x, gridx.y)
    Xfy, Yfy = numpy.meshgrid(gridy.x, gridy.y)

    for particle in particles:
        levelset_x = 0.5 - numpy.sqrt((Xfx-particle.x[0,0])**2 + (Yfx-particle.x[0,1])**2)
        levelset_y = 0.5 - numpy.sqrt((Xfy-particle.x[0,0])**2 + (Yfy-particle.x[0,1])**2)

    gridx.set_values(ibmf, levelset_x.transpose())
    gridy.set_values(ibmf, levelset_y.transpose()) 

    return
