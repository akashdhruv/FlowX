import numpy

def force_flow(gridc, gridx, gridy, scalars, particles, ibmf, ibmx, ibmy, velc, options):

    """
    Subroutine to compute forces on the fluid due to the presence of the immersed boundary
 
    Arguments
    ---------
    gridc : object
      Grid object for center variables

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

    u = gridx[velc][0,0,:,:]
    v = gridy[velc][0,0,:,:]
    
    ibc = gridc[ibmf][0,0,:,:]

    ibx = (ibc[:,:-1]+ibc[:,1:])/2.0
    iby = (ibc[:-1,:]+ibc[1:,:])/2.0

    indx = numpy.where(ibx >= 0.0)
    indy = numpy.where(iby >= 0.0)

    for particle in particles:
        u[indx] = particle.vel[0,0]
        v[indy] = particle.vel[0,1]

    return
