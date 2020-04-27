import numpy
from numba import jit

def advect_dynamic_grid(lmda,s,u,v,dt,dx,dy,nx,ny):
    """
    Subroutine to advect dynamic x-y grid

    """

    pfl = lmda>=0.0
    advectWENO3(s,u,v,dt,dx,dy,nx,ny)
    s = pfl*s

    return

def advect_solid(s,u,v,dt,dx,dy,nx,ny):
    """
    Subroutine to advect solid interface

    """
    advectWENO3(s,u,v,dt,dx,dy,nx,ny)

    return

def advectWENO3(s,u,v,dt,dx,dy,nx,ny):
    """
    Subroutine to add additional guard cells for WENO3 stencil

    """

    nguard_add = 2

    s_weno = numpy.zeros((nx+2*nguard_add,ny+2*nguard_add), dtype=float)
    u_weno = numpy.zeros((nx+2*nguard_add-1,ny+2*nguard_add), dtype=float)
    v_weno = numpy.zeros((nx+2*nguard_add,ny+2*nguard_add-1), dtype=float)

    s_weno[2:-2,2:-2] = s
    u_weno[2:-2,2:-2] = u
    v_weno[2:-2,2:-2] = v

    # x low
    s_weno[1,:] = 2*s_weno[2,:] - s_weno[3,:]
    s_weno[0,:] = 2*s_weno[1,:] - s_weno[2,:]

    # x high
    s_weno[-2,:] = 2*s_weno[-3,:] - s_weno[-4,:]
    s_weno[-1,:] = 2*s_weno[-2,:] - s_weno[-3,:]

    # y low
    s_weno[:,1] = 2*s_weno[:,2] - s_weno[:,3]
    s_weno[:,0] = 2*s_weno[:,1] - s_weno[:,2]

    # y high
    s_weno[:,-2] = 2*s_weno[:,-3] - s_weno[:,-4]
    s_weno[:,-1] = 2*s_weno[:,-2] - s_weno[:,-3]

    jit_advectWENO3(s_weno,u_weno,v_weno,dt,dx,dy,nx+2*nguard_add,ny+2*nguard_add)

    s[:,:] = s_weno[2:-2,2:-2]

    return

@jit(nopython=True)
def jit_advectWENO3(s,u,v,dt,dx,dy,nx,ny):

    so = numpy.copy(s)
    eps = 1E-15
   
    for i in range(3,nx-3):
        for j in range(3,ny-3):

            #--Velocities on faces used for divergence --> div(u*phi)

            ul = u[i-1,j]
            ur = u[i,j]
            vl = v[i,j-1]
            vr = v[i,j]

            #- WENO3 stencil in X direction --------------

            if (ur > 0):    # u = (+) Downwind

                s1r = so[i-2,j]
                s2r = so[i-1,j]
                s3r = so[i,j]
                s4r = so[i+1,j]
                s5r = so[i+2,j]

                rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. \
                      +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
                rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. \
                      +  1./4. *(    s2r           -    s4r )**2.
                rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. \
                      +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

                aT1r = 1./10. /  ( eps + rIS1r )**2.
                aT2r = 6./10. /  ( eps + rIS2r )**2.
                aT3r = 3./10. /  ( eps + rIS3r )**2.

                a1r = aT1r / ( aT1r + aT2r +aT3r )
                a2r = aT2r / ( aT1r + aT2r +aT3r )
                a3r = aT3r / ( aT1r + aT2r +aT3r )

                fT1r =  2./6.*s1r - 7./6.*s2r + 11./6.*s3r
                fT2r = -1./6.*s2r + 5./6.*s3r +  2./6.*s4r
                fT3r =  2./6.*s3r + 5./6.*s4r -  1./6.*s5r

            else:             # u = (-) Upwind
                     
                s1r = so[i-1,j]
                s2r = so[i,j]
                s3r = so[i+1,j]
                s4r = so[i+2,j]
                s5r = so[i+3,j]
  
                rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. \
                      +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
                rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. \
                      +  1./4. *(    s2r           -    s4r )**2.
                rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. \
                      +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

                aT1r = 3./10. /  ( eps + rIS1r )**2.
                aT2r = 6./10. /  ( eps + rIS2r )**2.
                aT3r = 1./10. /  ( eps + rIS3r )**2.

                a1r = aT1r / ( aT1r + aT2r +aT3r )
                a2r = aT2r / ( aT1r + aT2r +aT3r )
                a3r = aT3r / ( aT1r + aT2r +aT3r )

                fT1r = -1./6.*s1r + 5./6.*s2r +  2./6.*s3r
                fT2r =  2./6.*s2r + 5./6.*s3r -  1./6.*s4r
                fT3r =  11./6.*s3r - 7./6.*s4r + 2./6.*s5r

            if (ul > 0):       # u = (+) Downwind  

                s1l = so[i-3,j]
                s2l = so[i-2,j]
                s3l = so[i-1,j]
                s4l = so[i,j]
                s5l = so[i+1,j]
 
                rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. \
                      +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
                rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. \
                      +  1./4. *(    s2l           -    s4l )**2.
                rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. \
                      +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

                aT1l = 1./10. /  ( eps + rIS1l )**2.
                aT2l = 6./10. /  ( eps + rIS2l )**2.
                aT3l = 3./10. /  ( eps + rIS3l )**2.

                a1l = aT1l / ( aT1l + aT2l +aT3l )
                a2l = aT2l / ( aT1l + aT2l +aT3l )
                a3l = aT3l / ( aT1l + aT2l +aT3l )

                fT1l =  2./6.*s1l - 7./6.*s2l + 11./6.*s3l
                fT2l = -1./6.*s2l + 5./6.*s3l +  2./6.*s4l
                fT3l =  2./6.*s3l + 5./6.*s4l -  1./6.*s5l

            else:                # u = (-) Upwind

                s1l = so[i-2,j]
                s2l = so[i-1,j]
                s3l = so[i,j]
                s4l = so[i+1,j]
                s5l = so[i+2,j]   

                rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. \
                      +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
                rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. \
                      +  1./4. *(    s2l           -    s4l )**2.
                rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. \
                      +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

                aT1l = 3./10. /  ( eps + rIS1l )**2.
                aT2l = 6./10. /  ( eps + rIS2l )**2.
                aT3l = 1./10. /  ( eps + rIS3l )**2.

                a1l = aT1l / ( aT1l + aT2l +aT3l )
                a2l = aT2l / ( aT1l + aT2l +aT3l )
                a3l = aT3l / ( aT1l + aT2l +aT3l )

                fT1l = -1./6.*s1l + 5./6.*s2l +  2./6.*s3l
                fT2l =  2./6.*s2l + 5./6.*s3l -  1./6.*s4l
                fT3l =  11./6.*s3l - 7./6.*s4l + 2./6.*s5l

            # WENO3 interpolated PHI values at cell face... 
            frx = a1r*fT1r + a2r*fT2r + a3r*fT3r
            flx = a1l*fT1l + a2l*fT2l + a3l*fT3l

 
            #- WENO3 stencil in Y direction --------------

            if (vr > 0):     # u = (+) Downwind

                s1r = so[i,j-2]
                s2r = so[i,j-1]
                s3r = so[i,j]
                s4r = so[i,j+1]
                s5r = so[i,j+2]

                rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. \
                      +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
                rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. \
                      +  1./4. *(    s2r           -    s4r )**2.
                rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. \
                      +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

                aT1r = 1./10. /  ( eps + rIS1r )**2.
                aT2r = 6./10. /  ( eps + rIS2r )**2.
                aT3r = 3./10. /  ( eps + rIS3r )**2.

                a1r = aT1r / ( aT1r + aT2r +aT3r )
                a2r = aT2r / ( aT1r + aT2r +aT3r )
                a3r = aT3r / ( aT1r + aT2r +aT3r )

                fT1r =  2./6.*s1r - 7./6.*s2r + 11./6.*s3r
                fT2r = -1./6.*s2r + 5./6.*s3r +  2./6.*s4r
                fT3r =  2./6.*s3r + 5./6.*s4r -  1./6.*s5r

            else:                      # u = (-) Upwind

                s1r = so[i,j-1]
                s2r = so[i,j]
                s3r = so[i,j+1]
                s4r = so[i,j+2]
                s5r = so[i,j+3]

                rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. \
                      +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
                rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. \
                      +  1./4. *(    s2r           -    s4r )**2.
                rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. \
                      +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

                aT1r = 3./10. /  ( eps + rIS1r )**2.
                aT2r = 6./10. /  ( eps + rIS2r )**2.
                aT3r = 1./10. /  ( eps + rIS3r )**2.

                a1r = aT1r / ( aT1r + aT2r +aT3r )
                a2r = aT2r / ( aT1r + aT2r +aT3r )
                a3r = aT3r / ( aT1r + aT2r +aT3r )

                fT1r = -1./6.*s1r + 5./6.*s2r +  2./6.*s3r
                fT2r =  2./6.*s2r + 5./6.*s3r -  1./6.*s4r
                fT3r =  11./6.*s3r - 7./6.*s4r + 2./6.*s5r

            if (vl > 0):     # u = (+) Downwind

                s1l = so[i,j-3]
                s2l = so[i,j-2]
                s3l = so[i,j-1]
                s4l = so[i,j]
                s5l = so[i,j+1]

                rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. \
                      +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
                rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. \
                      +  1./4. *(    s2l           -    s4l )**2.
                rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. \
                      +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

                aT1l = 1./10. /  ( eps + rIS1l )**2.
                aT2l = 6./10. /  ( eps + rIS2l )**2.
                aT3l = 3./10. /  ( eps + rIS3l )**2.

                a1l = aT1l / ( aT1l + aT2l +aT3l )
                a2l = aT2l / ( aT1l + aT2l +aT3l )
                a3l = aT3l / ( aT1l + aT2l +aT3l )

                fT1l =  2./6.*s1l - 7./6.*s2l + 11./6.*s3l
                fT2l = -1./6.*s2l + 5./6.*s3l +  2./6.*s4l
                fT3l =  2./6.*s3l + 5./6.*s4l -  1./6.*s5l

            else:              # u = (-) Upwind

                 s1l = so[i,j-2]
                 s2l = so[i,j-1]
                 s3l = so[i,j]
                 s4l = so[i,j+1]
                 s5l = so[i,j+2]

                 rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. \
                       +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
                 rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. \
                       +  1./4. *(    s2l           -    s4l )**2.
                 rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. \
                       +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

                 aT1l = 3./10. /  ( eps + rIS1l )**2.
                 aT2l = 6./10. /  ( eps + rIS2l )**2.
                 aT3l = 1./10. /  ( eps + rIS3l )**2.

                 a1l = aT1l / ( aT1l + aT2l +aT3l )
                 a2l = aT2l / ( aT1l + aT2l +aT3l )
                 a3l = aT3l / ( aT1l + aT2l +aT3l )

                 fT1l = -1./6.*s1l + 5./6.*s2l +  2./6.*s3l
                 fT2l =  2./6.*s2l + 5./6.*s3l -  1./6.*s4l
                 fT3l =  11./6.*s3l - 7./6.*s4l + 2./6.*s5l

            #  WENO3 interpolated PHI values at cell face... 
            fry = a1r*fT1r + a2r*fT2r + a3r*fT3r
            fly = a1l*fT1l + a2l*fT2l + a3l*fT3l
                
            s[i,j] = so[i,j] - dt*(frx*ur - flx*ul)/dx - dt*(fry*vr - fly*vl)/dy

    return
