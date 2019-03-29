import numpy

def divergence(gridc, gridx, gridy, ivar, dvar):
    """Solve the divergence 
    
    Arguments
    ---------
    gridc : grid object
        Grid containing data
    gridx : 
    
    gridy :
        
    ivar : string
        Name of the grid variable of the u and v data
    dvar : string
        Name of the grid variable of the divergence
    
    """
    
    u_ivar = gridx.get_variable_indices(ivar)

    v_ivar = gridy.get_variable_indices(ivar)
    
    div_ivar = gridc.get_variable_indices(dvar) 
    dx, dy = gridc.dx, gridc.dy
    
    Udata = gridx.data[:,:,u_ivar]
    Vdata = gridy.data[:,:,v_ivar]

    gridc.data[1:-1, 1:-1, div_ivar] = ((Udata[1:, 1:-1] - Udata[:-1, 1:-1])/dx 
                                     +  (Vdata[1:-1, 1:] - Vdata[1:-1, :-1])/dy)
  
    gridc.fill_guard_cells(dvar)

    return
