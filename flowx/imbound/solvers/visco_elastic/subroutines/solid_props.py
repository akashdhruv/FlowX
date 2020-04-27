import numpy
from scipy import special

def solid_props(s,xmus,mu_solid,dx,dy,nx,ny):

    eps = 1.E-10
       
    psi  = (1. + special.erf(-s/(2.*dx)))/2.

    xmus[:,:] = psi*(mu_solid-0.) + 0.

    return
