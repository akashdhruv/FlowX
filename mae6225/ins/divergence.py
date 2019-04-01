import numpy

def divergence(gridc, gridx, gridy, ivar, dvar, ifac=1.):
    """Divergence operator
    
    Arguments
    ---------
    gridc : grid object (center)
        Grid containing data in on cell centers

    gridx : grid object (x-direction)
        Grid containing data in x-direction

    gridy : grid object (y-direction)
        Grid containing data in y-direction
       
    ivar : string
        Name of the face-centered grid variable for velocity

    dvar : string
        Name of the cell-centered grid variable to store divergence

    ifac : float (optional)
	 Multiplying factor for time-step
    
    """

    u   = gridx.get_values(ivar)
    v   = gridy.get_values(ivar)

    div = gridc.get_values(dvar)
    
    dx, dy = gridc.dx, gridc.dy
    

    div[1:-1,1:-1] = ((u[1:, 1:-1] - u[:-1, 1:-1])/dx
                     +(v[1:-1, 1:] - v[1:-1, :-1])/dy)/ifac

  
    gridc.fill_guard_cells(dvar)

    return
