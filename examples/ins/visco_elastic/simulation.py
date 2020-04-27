"""User defined module for simulation."""
  
import numpy


def get_initial(gridc, gridx, gridy, ivar, pvar):
    """Compute and set the analytical solution.

    Arguments
    ---------
    gridc : flowx.Grid object
        Grid containing cell-center data.
    
    gridx : flowx.Grid object
        Grid containing x-face data

    girdy : flowx.Grid object
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

    u[:,:] =  -0.1*numpy.pi*numpy.sin(2*numpy.pi*(Xfx+0.5))*numpy.cos(2*numpy.pi*(Yfx+0.5))
    v[:,:] =   0.1*numpy.pi*numpy.sin(2*numpy.pi*(Yfy+0.5))*numpy.cos(2*numpy.pi*(Xfy+0.5))

    return
