import numpy


def diffusion(grid, ivar, alpha):
    """Diffusion operator
    
    Arguments
    ---------
    grid  : grid object
        Grid containing data for a given stencil

    ivar  : string
        Name of the grid variable to be operated on

    alpha : float
          Diffusion coefficient

    Returns
    --------
    D : float
      Array contains value of diffusion terms
          
    """

    f = grid.get_values(ivar)

    dx, dy = grid.dx, grid.dy

    D = alpha*((f[2:,1:-1] - 2.0*f[1:-1,1:-1] + f[:-2,1:-1])/(dx**2) 
              +(f[1:-1,2:] - 2.0*f[1:-1,1:-1] + f[1:-1,:-2])/(dy**2))
                    
    return D
