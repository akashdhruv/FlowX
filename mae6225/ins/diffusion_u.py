import numpy


def diffusion_u(gridx, ivar, dvar, Re):
    nu = 1/Re
    i_ivar, i_dvar = gridx.get_variable_indices(ivar, dvar)
    dx, dy = gridx.dx, gridx.dy
    gridx.data[1:-1, 1:-1, i_dvar] = nu*((gridx.data[2:,1:-1,i_ivar]-2.0*gridx.data[1:-1,1:-1,i_ivar] 
                                     + gridx.data[:-2,1:-1,i_ivar])/(dx**2) 
                                     + (gridx.data[1:-1,2:,i_ivar] - 2.0*gridx.data[1:-1,1:-1,i_ivar] 
                                     + gridx.data[1:-1,:-2,i_ivar])/(dy**2))
                    
return