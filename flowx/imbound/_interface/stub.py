def force_flow(gridc, gridx, gridy, scalars, particles, ibmf, ibmx, ibmy, velc, options):
    """
    Stub subroutine to compute forces
 
    Arguments
    ---------
    gridc : object
      Grid object for cell center variables

    gridx : object
      Grid object for x-face variables

    gridy : object
      Grid object for y-face variables

    scalars: object
       Scalars object to access time-step and Reynold number

    particles: object
       Particles object to access time-step 

    ibmf : string for forcing variable

    velc : string for velocity variable
    """
    pass

def map_to_grid(gridc, particles, ibmf, ibmx, ibmy, search_function, options):

    """
    Stub subroutine to compute IB mapping on grid
 
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

    return None

def advect(gridc, gridx, gridy, scalars, ibmf, ibmx, ibmy, velc, options):

    """
    Stub subroutine for level set advection
    """

    return
