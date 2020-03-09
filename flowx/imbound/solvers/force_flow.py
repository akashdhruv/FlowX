import numpy

def force_flow_stub(gridx, gridy, scalars, particles, ibmf, velc):

    """
    Stub subroutine to compute forces
 
    Arguments
    ---------
    gridx : object
      Grid object for x-face variables

    gridy : object
      Grid object for y-face variables

    scalars: object
       Scalars object to access time-step and Reynold number

    particles: object
       Particles object to access time-step 

    ibmf : string for forcing variable

    velc : string for velocity variable
    """

    return
 
def force_flow_levelset(gridx, gridy, scalars, particles, ibmf, velc):

    """
    Subroutine to compute forces on the fluid due to the presence of the immersed boundary
 
    Arguments
    ---------
    gridx : object
      Grid object for x-face variables

    gridy : object
      Grid object for y-face variables

    scalars: object
       Scalars object to access time-step and Reynold number

    particles: object
       Object containing immersed boundary information

    ibmf : string for forcing variable

    velc : string for velocity variable
    """

    Xfx, Yfx = numpy.meshgrid(gridx.x, gridx.y)
    Xfy, Yfy = numpy.meshgrid(gridy.x, gridy.y)

    Xfx = Xfx.transpose()
    Yfx = Yfx.transpose()
    
    Xfy = Xfy.transpose()
    Yfy = Yfy.transpose()

    u = gridx.get_values(velc)
    v = gridy.get_values(velc)

    ibx = gridx.get_values(ibmf)
    iby = gridy.get_values(ibmf)

    for particle in particles:
        indx = numpy.where(ibx >= 0.0)
        indy = numpy.where(iby >= 0.0)
        u[indx] = particle.vel[0,0]
        v[indy] = particle.vel[0,1]
    return
