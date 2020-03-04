"""Implementation of the immersed boundary module"""

from flowx.imbound.imbound_interface import imbound_interface

class imbound_main(imbound_interface):

    def __init__(self, imbound_vars=None, **kwargs):

        """
        Constructor for the imbound unit

        Arguments
        ---------

        imbound_vars : list
                List of string for field variables required by ins unit
               
                imbound_vars[0] --> indicator variable for immersed boundary
                imbound_vars[1] --> velocity variable

        **kwargs : Dictionary of keyword arguments

        'with_ib' - keyword to indicate if immersed boundary is present or not
                  - True if present
                  - False if not present

        kwargs['with_ib'] = False --> default
                          = True

        Error generated if kwargs['with_ib'] is True and imbound_vars is None

        """

        from flowx.imbound.solvers.force_flow import force_flow_stub, force_flow_levelset
        from flowx.imbound.solvers.map_to_grid import map_to_grid_stub, map_to_grid_levelset

        self._ibmf = 'stub'
        self._velc = 'stub'

        if imbound_vars is not None:
            self._ibmf = imbound_vars[0]
            self._velc = imbound_vars[1] 


        self._with_ib = False

        if 'with_ib' in kwargs: self._with_ib = kwargs['with_ib']

        self._force_flow  = force_flow_stub
        self._map_to_grid = map_to_grid_stub
 
        if self._with_ib:
            self._force_flow = force_flow_levelset
            self._map_to_grid = map_to_grid_levelset

        if(self._with_ib and imbound_vars is None): raise ValueError('imbound_vars cannot be None when with_ib is True')

        return
 
    def map_to_grid(self, gridx, gridy, particles):
        """
        Subroutine to map immersed boundary on grid
 
        Arguments
        ---------
        gridx : object
          Grid object for x-face variables

        gridy : object
          Grid object for y-face variables

        particles: object
           Object containing immersed boundary information
        """

        self._map_to_grid(gridx, gridy, particles, self._ibmf, self._velc)

        return

    def force_flow(self, gridx, gridy, scalars, particles):

        """
        Subroutine to compute immersed boundary forces
 
        Arguments
        ---------
        gridx : object
          Grid object for x-face variables

        gridy : object
          Grid object for y-face variables

        scalars: object
           Scalars object to access time-step and Reynold number

        particles: object
           Object containing immersed boundary information
        """

        self._force_flow(gridx, gridy, scalars, particles, self._ibmf, self._velc)

        return
