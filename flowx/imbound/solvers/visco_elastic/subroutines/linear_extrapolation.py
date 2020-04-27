import numpy
from numba import jit

def linear_extrapolation(lmda,s,sn,u,v,dx,dy,nx,ny):

    pfl = lmda<0.0
    so = numpy.copy(s)

    delta_t = dx/2.
    
    up = u[1:-1,1:-1]
    vp = v[1:-1,1:-1]

    sxplus = (so[2:,1:-1]-so[1:-1,1:-1])/dx
    sxmins = (so[1:-1,1:-1]-so[:-2,1:-1])/dx

    syplus = (so[1:-1,2:]-so[1:-1,1:-1])/dy
    symins = (so[1:-1,1:-1]-so[1:-1,:-2])/dy

    jit_linear_extrapolation(pfl,s,so,sn,up,vp,sxplus,sxmins,syplus,symins,delta_t)

    return

@jit(nopython=True)
def jit_linear_extrapolation(pfl,s,so,sn,up,vp,sxplus,sxmins,syplus,symins,dt):

    up_zeros = numpy.zeros_like(up)
    vp_zeros = numpy.zeros_like(vp)

    s[1:-1,1:-1] = so[1:-1,1:-1] + pfl[1:-1,1:-1]*dt*(sn[1:-1,1:-1] \
                                 - numpy.maximum(up,up_zeros)*sxmins - numpy.minimum(up,up_zeros)*sxplus \
                                 - numpy.maximum(vp,vp_zeros)*symins - numpy.minimum(vp,vp_zeros)*syplus )

    return
