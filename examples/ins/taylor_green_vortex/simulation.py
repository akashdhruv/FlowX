"""User defined module for simulation."""
  
import numpy


def get_initial(gridc, gridx, gridy, ivar, pvar):
    """Compute and set the analytical solution.

    Arguments
    ---------
    gridc : mae6225.Grid object
        Grid containing cell-center data.
    
    gridx : mae6225.Grid object
        Grid containing x-face data

    girdy : mae6225.Grid object
        Grid containing y-face data

    ivar : string
        Name of the velocity variable on the grid.

    pvar : string
        Name of the pressure variable on the grid.

    """

    Xcc, Ycc = numpy.meshgrid(gridc.x, gridc.y)
    Xfx, Yfx = numpy.meshgrid(gridx.x, gridx.y)
    Xfy, Yfy = numpy.meshgrid(gridy.x, gridy.y)


    Xcc = Xcc.transpose()
    Ycc = Ycc.transpose()

    Xfx = Xfx.transpose()
    Yfx = Yfx.transpose()

    Xfy = Xfy.transpose()
    Yfy = Yfy.transpose()
 
    u = gridx.get_values(ivar)
    v = gridy.get_values(ivar)
    p = gridc.get_values(pvar)

    u[:,:] =  -numpy.cos(Xfx)*numpy.sin(Yfx)
    v[:,:] =   numpy.sin(Xfy)*numpy.cos(Yfy)
    p[:,:] = -(numpy.cos(2*Xcc)+numpy.sin(2*Ycc))/4

    return

def get_analytical(gridc, gridx, gridy, asol, ifac):
    """Compute and set the analytical solution.

    Arguments
    ---------
    gridc : mae6225.Grid object
        Grid containing cell-center data.
    
    gridx : mae6225.Grid object
        Grid containing x-face data

    girdy : mae6225.Grid object
        Grid containing y-face data

    asol : string
        Name of the variable on the grid.

    ifac : float
       Time step

    """

    Xcc, Ycc = numpy.meshgrid(gridc.x, gridc.y)
    Xfx, Yfx = numpy.meshgrid(gridx.x, gridx.y)
    Xfy, Yfy = numpy.meshgrid(gridy.x, gridy.y)

    Xcc = Xcc.transpose()
    Ycc = Ycc.transpose()

    Xfx = Xfx.transpose()
    Yfx = Yfx.transpose()

    Xfy = Xfy.transpose()
    Yfy = Yfy.transpose()

    u = gridx.get_values(asol)
    v = gridy.get_values(asol)
    p = gridc.get_values(asol)

    u[:,:] =  -numpy.exp(-2*ifac)*numpy.cos(Xfx)*numpy.sin(Yfx)
    v[:,:] =   numpy.exp(-2*ifac)*numpy.sin(Xfy)*numpy.cos(Yfy)
    p[:,:] =  -numpy.exp(-4*ifac)*(numpy.cos(2*Xcc)+numpy.sin(2*Ycc))/4

    return
