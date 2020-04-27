import numpy
import time

from flowx.imbound.solvers.visco_elastic.subroutines import advect_dynamic_grid, advect_solid, \
                                                            redistance_solid, normal_vector_solid, \
                                                            directional_derivative, constant_extrapolation, \
                                                            linear_extrapolation

def advect_visco(gridc, gridx, gridy, scalars, ibmf, ibmx, ibmy, velc, options):

    """

    """

    nx, ny = gridc.nx, gridc.ny
    dx, dy = gridc.dx, gridc.dy
    dt = scalars.dt
    lset_iter = options['lset_iter']
    extp_iter = options['extp_iter']

    phi = gridc.get_values(ibmf)
    lmx = gridc.get_values(ibmx)
    lmy = gridc.get_values(ibmy)

    adfx = numpy.zeros_like(phi)
    adfy = numpy.zeros_like(phi)
    ddsn = numpy.zeros_like(phi)

    u = gridx.get_values(velc)
    v = gridy.get_values(velc)

    #--------------Advect dynamic X-Y grid---------------
    advect_dynamic_grid(phi,lmx,u,v,dt,dx,dy,nx+2,ny+2)
    advect_dynamic_grid(phi,lmy,u,v,dt,dx,dy,nx+2,ny+2)

    #---------Find normal vectors---------------------
    normal_vector_solid(phi,adfx,adfy,dx,dy,nx+2,ny+2)

    #---------Extrapolate X grid----------------------
    directional_derivative(phi,lmx,ddsn,adfx,adfy,dx,dy,nx+2,ny+2)
    
    for _iter in range(extp_iter): 
        constant_extrapolation(phi,ddsn,adfx,adfy,dx,dy,nx+2,ny+2)

    for _iter in range(extp_iter):
        linear_extrapolation(phi,lmx,ddsn,adfx,adfy,dx,dy,nx+2,ny+2)    
        gridc.fill_guard_cells(ibmx)

    #---------Extrapolate Y grid----------------------
    directional_derivative(phi,lmy,ddsn,adfx,adfy,dx,dy,nx+2,ny+2)

    for _iter in range(extp_iter): 
        constant_extrapolation(phi,ddsn,adfx,adfy,dx,dy,nx+2,ny+2)

    for _iter in range(extp_iter):
        linear_extrapolation(phi,lmy,ddsn,adfx,adfy,dx,dy,nx+2,ny+2)    
        gridc.fill_guard_cells(ibmy)

    #--------Advect solid interface-------------------
    advect_solid(phi,u,v,dt,dx,dy,nx+2,ny+2) 
    gridc.fill_guard_cells(ibmf)

    #--------Redistance solid interface---------------
    phi_old = numpy.copy(phi)
    lsDT = numpy.sqrt(dx**2 + dy**2)/2.

    for _iter in range(lset_iter): 
        redistance_solid(phi,phi_old,lsDT,dx,dy,nx+2,ny+2)
        gridc.fill_guard_cells(ibmf)

    return
