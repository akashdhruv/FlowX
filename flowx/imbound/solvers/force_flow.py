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

    velx = gridx.get_values(velc)
    vely = gridy.get_values(velc)

    ibx = gridx.get_values(ibmf)
    iby = gridy.get_values(ibmf)

    indx = numpy.where(ibx >= 0.0)
    indy = numpy.where(iby >= 0.0)

    velx[indx] = particles[0].velx
    vely[indy] = particles[0].vely

    return
