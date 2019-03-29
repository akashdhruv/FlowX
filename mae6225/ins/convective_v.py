import numpy

def convective_v(gridx, gridy, ivar, cvar)
    """ Solving convection in the y-direction for incompressible NS Equation.
    
    Arguments
    ---------
    gridx : grid object (x-direction)
        Grid containing data in x-direction
    gridy: grid object (ydirection)
        Grid containing data in y-direction
    ivar: string
        Name of the grid variable of the numerical solution
    cvar: string
        Name of the grid variable of the convection term
    
        
    Returns
    -------
    ites : integer
        Number of iterations computed
    residual : float
        Final residual
    verbose : bool, optional
        Set True to display convergence information;
        Default : False
        
    """
    
    i_ivar, i_cvar= grid.get_variable_indices(ivar, cvar)
    dx, dy = gridy.dx, gridy.dy
    nx, ny = gridy.nx, gridy.ny
    C = numpy.zeros((nx+2,ny+1))
    # Read u-velocity value from the grid object
    for i in range (1,nx+1):
        for j in range (1, ny):
            uij = gridx.data[i, j, i_var]
      
            uij1 = gridx.data[i, j+1, i_var]
    
            ui_1j = gridx.data[i-1, j, i_var]
      
            ui_1j1 = gridx.data[i-1, j+1, i_var]
    
     # Read v-velocity value from the grid object
            vij = gridy.data[i, j, i_var]
    
            vi1j = gridy.data[i+1, j, i_var]
    
            vij1 = gridy.data[i, j+1, i_var]
    
            vi_1j = gridy.data[i-1, j, i_var]
    
            vij_1 = gridy.data[i, j-1, i_var]
    
    #Read u*v-velocity value from the northeast-corner and northwest-corner on grid object
    
            uv_ne = ((uij + uij1)/2)*(((vij + vi1j)/2))
    
            uv_nw = ((ui_1j + ui_1j1)/2)*(((vij + vi_1j)/2))
    
    #uv_bot = [(uij + uij_1)/2]*[((vij + vi1j)/2)]
    
    # Compute u*v-velocity value with respect to x grid space
    
            duv_dx = ((uv_nw - uv_ne)/dx)
    
    #Compute v-velocity value squared  with respect to y grid space
    
            v_n = ((vij+vij1)/2)
            v_s = ((vij+vij_1)/2)
    
            dv_2_dy = (((v_n**2)-(v_s**2))/dy)
    
    #Compute convection-v value in y-direction
    
            C[i,j] = ((-dv_2_dx) - (duv_dx))
    
    # store the convection-term values in the grid object
    
    gridy.set_values(cvar,C)
    
    # return nothing as we store directly the data in the grid object
    return
