import numpy
import time

def map_to_grid(gridc, particles, ibmf, ibmx, ibmy, search_function, options):

    """
    Subroutine to compute IB mapping on grid using the level set function
 
    Arguments
    ---------
    gridc : object
      Grid object for x-face variables

    particles: object
       Object containing immersed boundary information

    ibmf : string for forcing variable

    search_function : 
         Search function

    """
    for particle in particles:
        ites = search_function(gridc,particle,ibmf,options)

    if (options['verbose']):
        print("Mapping Iterations: ", ites)

    return ites
