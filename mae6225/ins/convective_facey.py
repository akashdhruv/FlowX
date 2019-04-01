import numpy

def convective_facey(gridx, gridy, ivar):
    """Convection operator for the facey grid
    
    Arguments
    ---------
    gridx : grid object (x-direction)
        Grid containing data in x-direction

    gridy : grid object (y-direction)
        Grid containing data in y-direction

    ivar: string
        Name of the grid variable to be operated on

    Returns
    --------
    C : float 
      Array contains value of convective terms
          
    """
    
    fx = gridx.get_values(ivar)
    fy = gridy.get_values(ivar)

    dx = gridx.dx
    dy = gridy.dy

    nx, ny = gridy.nx, gridy.ny

    C = numpy.zeros((nx+2,ny+1))

    for i in range (1,nx+1):
        for j in range (1, ny):

         # Get stencil in x-dir
            fxij    = fx[i, j]
      
            fxij1   = fx[i, j+1]
    
            fxi_1j  = fx[i-1, j]
      
            fxi_1j1 = fx[i-1, j+1]
    
        # Get stencil in y-dir
            fyij    = fy[i, j]
    
            fyi1j   = fy[i+1, j]
    
            fyij1   = fy[i, j+1]
    
            fyi_1j  = fy[i-1, j]
    
            fyij_1  = fy[i, j-1]
    
        #Read fx*fy value from the northeast-corner and northwest-corner on grid object
    
            fxfy_ne = ((fxij   + fxij1)/2)   * ((fyij + fyi1j)/2)
    
            fxfy_nw = ((fxi_1j + fxi_1j1)/2) * ((fyij + fyi_1j)/2)

        #Read fy*fy value from the nort and south on grid object
   
            fyfy_n  = ((fyij+fyij1)/2)**2
   
            fyfy_s  = ((fyij+fyij_1)/2)**2
   
        # Compute derivates
    
            dfxfy_dx = (fxfy_nw - fxfy_ne)/dx   

            dfyfy_dy = (fyfy_n  - fyfy_s)/dy
    
        #Compute convection-v value in y-direction
    
            C[i,j] = - dfxfy_dx - dfyfy_dy
       
   
    # return only internal values for C
    return C[1:-1,1:-1]
