import numpy
from numba import jit

def redistance_solid(s,soo,dt,dx,dy,nx,ny):

    _jit_redistance_solid(s,soo,dt,dx,dy,nx,ny)

    return

@jit(nopython=True)
def _jit_redistance_solid(s,soo,dt,dx,dy,nx,ny):

    eps = 1E-14
        
    sgn = soo/numpy.abs(soo+eps)
    so  = numpy.copy(s)

    s[:,:] = 0.

    for i in range(1,nx-1):
        for j in range(1,ny-1):

            sm  = so[i,j]

            sxl = so[i-1,j]
            sxr = so[i+1,j]
            syl = so[i,j-1]
            syr = so[i,j+1]

            if ((soo[i,j]*soo[i-1,j] < 0) or \
                (soo[i,j]*soo[i+1,j] < 0) or \
                (soo[i,j]*soo[i,j-1] < 0) or \
                (soo[i,j]*soo[i,j+1] < 0)):

                s[i,j] = soo[i,j]

            else:

                #- First Order Upwind -----
                ap = numpy.maximum((sm-sxl),0.)/dx
                an = numpy.minimum((sm-sxl),0.)/dx

                bp = numpy.maximum((sxr-sm),0.)/dx
                bn = numpy.minimum((sxr-sm),0.)/dx

                cp = numpy.maximum((sm-syl),0.)/dy
                cn = numpy.minimum((sm-syl),0.)/dy

                dp = numpy.maximum((syr-sm),0.)/dy
                dn = numpy.minimum((syr-sm),0.)/dy

                if (so[i,j] > 0.):
                    agf = ( numpy.sqrt(numpy.maximum(ap**2,bn**2) + numpy.maximum(cp**2,dn**2)) ) -1.0

                elif (so[i,j] < 0.):
                    agf = ( numpy.sqrt(numpy.maximum(an**2,bp**2) + numpy.maximum(cn**2,dp**2)) ) -1.0

                else:
                    agf = 0.

                # - Solve Level Set Re-distance equation ---------

                s[i,j] = so[i,j] - dt*(sgn[i,j]*(agf))

    if(s[i,j]*soo[i,j] < 0.0): print("WARNING: LS Dist Function Changed Signs - ",i,j)

    return
