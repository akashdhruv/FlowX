import numpy
from numba import jit

def directional_derivative(lmda,so,ddsn,adfx,adfy,dx,dy,nx,ny):

    pfl = lmda>=0.0

    up = adfx[1:-1,1:-1]
    vp = adfy[1:-1,1:-1]

    up_zeros = numpy.zeros_like(up)
    vp_zeros = numpy.zeros_like(vp)

    sxplus = (so[2:,1:-1]-so[1:-1,1:-1])/dx
    sxmins = (so[1:-1,1:-1]-so[:-2,1:-1])/dx
 
    syplus = (so[1:-1,2:]-so[1:-1,1:-1])/dy
    symins = (so[1:-1,1:-1]-so[1:-1,:-2])/dy

    ddsn[1:-1,1:-1] = numpy.maximum(up,up_zeros)*sxmins + numpy.minimum(up,up_zeros)*sxplus + \
                      numpy.maximum(vp,vp_zeros)*symins + numpy.minimum(vp,vp_zeros)*syplus
       
    ddsn[0,:] = ddsn[1,:]
    ddsn[-1,:] = ddsn[-2,:]

    ddsn[:,0] = ddsn[:,1]
    ddsn[:,-1] = ddsn[:,-2]

    return
