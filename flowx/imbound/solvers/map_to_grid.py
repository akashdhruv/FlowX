import numpy
import time

def map_to_grid_stub(gridc, gridx, gridy, particles, ibmf, search_function, options):

    """
    Stub subroutine to compute IB mapping on grid
 
    Arguments
    ---------
    gridc : object
      Grid object for x-face variables

    gridx : object
      Grid object for x-face variables

    gridy : object
      Grid object for y-face variables

    particles: object
       Object containing immersed boundary information

    ibmf : string for forcing variable

    search_function :
           Search function
    """

    return None

def map_to_grid_levelset(gridc, gridx, gridy, particles, ibmf, search_function, options):

    """
    Subroutine to compute IB mapping on grid using the level set function
 
    Arguments
    ---------
    gridc : object
      Grid object for x-face variables

    gridx : object
      Grid object for x-face variables

    gridy : object
      Grid object for y-face variables

    particles: object
       Object containing immersed boundary information

    ibmf : string for forcing variable

    search_function : 
         Search function

    """

    X, Y = numpy.meshgrid(gridc.x, gridc.y)

    nx, ny = gridc.nx, gridc.ny
    dx, dy = gridc.dx, gridc.dy

    X = X.transpose()
    Y = Y.transpose()

    IBx = gridx.get_values(ibmf)
    IBy = gridy.get_values(ibmf)
    IBc = gridc.get_values(ibmf)

    for particle in particles:
        ites = search_function(X, Y, IBc, nx+2, ny+2, particle, options)

    IBx[:,:] = (IBc[:-1,:]+IBc[1:,:])/2.0
    IBy[:,:] = (IBc[:,:-1]+IBc[:,1:])/2.0

    if (options['verbose']):
        print("Mapping Iterations: ", ites)

    return ites
