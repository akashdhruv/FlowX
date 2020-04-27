import numpy

from flowx.imbound.solvers.visco_elastic.subroutines import solid_props, solid_stress, solid_ustar

def force_flow_visco(gridc, gridx, gridy, scalars, particles, ibmf, ibmx, ibmy, velc):

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

    nx, ny = gridc.nx, gridc.ny
    dx, dy = gridc.dx, gridc.dy
    dt = scalars.dt
    Re_s, mu_s = scalars.Re_s, scalars.mu_s

    u = gridx.get_values(velc)
    v = gridy.get_values(velc)
    
    phi = gridc.get_values(ibmf)
    lmx = gridc.get_values(ibmx)
    lmy = gridc.get_values(ibmy)

    xmus = numpy.zeros_like(phi)
    lms1 = numpy.zeros_like(phi)
    lms2 = numpy.zeros_like(phi)
    lms3 = numpy.zeros_like(phi)
    lms4 = numpy.zeros_like(phi)

    solid_props(phi,xmus,mu_s,dx,dy,nx+2,ny+2)

    solid_stress(phi,lmx,lmy,lms1,lms2,lms3,lms4,dx,dy,nx+2,ny+2)

    solid_ustar(u,v,xmus,lms1,lms2,lms3,lms4,Re_s,dt,dx,dy,nx+2,ny+2)

    return
