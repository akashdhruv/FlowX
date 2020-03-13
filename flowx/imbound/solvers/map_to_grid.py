import numpy
from numba import jit

def map_to_grid_stub(gridx, gridy, particles, ibmf):

    """
    Stub subroutine to compute IB mapping on grid
 
    Arguments
    ---------
    gridx : object
      Grid object for x-face variables

    gridy : object
      Grid object for y-face variables

    particles: object
       Object containing immersed boundary information

    ibmf : string for forcing variable

    """

    return

def map_to_grid_levelset(gridx, gridy, particles, ibmf):

    """
    Subroutine to compute IB mapping on grid using the level set function
 
    Arguments
    ---------
    gridx : object
      Grid object for x-face variables

    gridy : object
      Grid object for y-face variables

    particles: object
       Object containing immersed boundary information

    ibmf : string for forcing variable

    """

    Xfx, Yfx = numpy.meshgrid(gridx.x, gridx.y)
    Xfy, Yfy = numpy.meshgrid(gridy.x, gridy.y)

    nx, ny = gridx.nx, gridx.ny

    Xfx = Xfx.transpose()
    Yfx = Yfx.transpose()

    Xfy = Xfy.transpose()
    Yfy = Yfy.transpose()

    IBx = gridx.get_values(ibmf)
    IBy = gridy.get_values(ibmf)

    for particle in particles:
        map_classical_search(nx+1, ny+2, Xfx, Yfx, IBx, particle.x[1:,:], particle.nnp-1)
        map_classical_search(nx+2, ny+1, Xfy, Yfy, IBy, particle.x[1:,:], particle.nnp-1)

    gridx.set_values(ibmf, levelset_x.transpose())
    gridy.set_values(ibmf, levelset_y.transpose()) 

    return

@jit(nopython=True)
def map_classical_search(nx, ny, x, y, phi, points, np):

    for i in range(nx):
        for j in range(ny):
            for p in range(np):


    return

def map_grovers_search(nx, ny, x, y, phi, points, np):

    return
