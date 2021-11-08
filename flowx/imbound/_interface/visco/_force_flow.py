import numpy

from . import solid_props, solid_stress, solid_ustar, normal_vector_solid, constant_extrapolation

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

    nx, ny = gridc.nx, gridc.ny
    dx, dy = gridc.dx, gridc.dy
    dt = scalars.dt
    Re_s, mu_s = scalars.Re_s, scalars.mu_s

    extrap_iter = options['extrap_solid']

    u = gridx[velc]
    v = gridy[velc]
    
    phi = gridc[ibmf]
    lmx = gridc[ibmx]
    lmy = gridc[ibmy]

    xmus = numpy.zeros_like(phi)
    lms1 = numpy.zeros_like(phi)
    lms2 = numpy.zeros_like(phi)
    lms3 = numpy.zeros_like(phi)
    lms4 = numpy.zeros_like(phi)
    adfx = numpy.zeros_like(phi)
    adfy = numpy.zeros_like(phi)

    #----------Assign solid properties---------------
    solid_props(phi,xmus,mu_s,dx,dy,nx+2,ny+2)

    #----------Calculate solid stress terms----------
    solid_stress(phi,lmx,lmy,lms1,lms2,lms3,lms4,dx,dy,nx+2,ny+2)

    #---------Find normal vectors---------------------
    normal_vector_solid(phi,adfx,adfy,dx,dy,nx+2,ny+2)

    #---------Extrapolation of stress terms----------    
    for _iter in range(extrap_iter):
        constant_extrapolation(phi,lms1,adfx,adfy,dx,dy,nx+2,ny+2)
        constant_extrapolation(phi,lms2,adfx,adfy,dx,dy,nx+2,ny+2)
        constant_extrapolation(phi,lms3,adfx,adfy,dx,dy,nx+2,ny+2)
        constant_extrapolation(phi,lms4,adfx,adfy,dx,dy,nx+2,ny+2)
 
    #---------------Update velocity------------------
    solid_ustar(u,v,xmus,lms1,lms2,lms3,lms4,Re_s,dt,dx,dy,nx+2,ny+2)

    return
