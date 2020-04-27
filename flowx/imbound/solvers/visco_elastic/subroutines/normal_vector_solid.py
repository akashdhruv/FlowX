import numpy
from numba import jit

def normal_vector_solid(sd,adfx,adfy,dx,dy,nx,ny):

    sxl = sd[:-2,1:-1]
    sxr = sd[2:,1:-1]
    syl = sd[1:-1,:-2]
    syr = sd[1:-1,2:]

    adf = numpy.sqrt(((sxr-sxl)/2./dx)**2 + ((syr-syl)/2./dy)**2)

    adfx[1:-1,1:-1] = -(sxr-sxl)/2./dx / adf
    adfy[1:-1,1:-1] = -(syr-syl)/2./dy / adf

    adfx[0,:] = adfx[1,:]
    adfx[-1,:] = adfx[-2,:]

    adfy[:,0] = adfy[:,1]
    adfy[:,-1] = adfy[:,-2]

    return
