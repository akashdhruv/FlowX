"""Implementation of the immersed boundary module"""

from flowx.imbound.imbound_interface import imbound_interface

class imbound_main(imbound_interface):

    def __init__(self, imbound_vars=None, imbound_info=None):

        """
        Constructor for the imbound unit

        Arguments
        ---------

        imbound_vars : list
                List of string for field variables required by ins unit
               
                imbound_vars[0] --> indicator variable for immersed boundary
                imbound_vars[1] --> velocity variable

        imbound_info : Dictionary of keyword arguments

        'with_ib' - keyword to indicate if immersed boundary is present or not
                  - True if present
                  - False if not present

        imbound_info['with_ib'] = False --> default
                                = True

        Error generated if kwargs['with_ib'] is True and imbound_vars is empty

        """

        from flowx.imbound.solvers.force_flow import force_flow_stub, force_flow_levelset
        from flowx.imbound.solvers.map_to_grid import map_to_grid_stub, map_to_grid_levelset

        self._ibmf = 'stub'
        self._velc = 'stub'
 
        self._with_ib = False
 
        self._force_flow  = force_flow_stub
        self._map_to_grid = map_to_grid_stub

        if imbound_info and 'with_ib' in imbound_info: self._with_ib = imbound_info['with_ib']

        if self._with_ib:
            self._force_flow = force_flow_levelset
            self._map_to_grid = map_to_grid_levelset

            if imbound_vars:
                self._ibmf = imbound_vars[0]
                self._velc = imbound_vars[1] 
            else:
                raise ValueError('imbound_vars cannot be empty when body flag is true')

        else:
            print('Warning: Immersed Boundary unit is a stub, no forcing will be applied.')

        return
 
    def map_to_grid(self, domain_data_struct):
        """
        Subroutine to map immersed boundary on grid
 
        Arguments
        ---------
        domain_data_struct : object list

        """

        _gridx = domain_data_struct[1]
        _gridy = domain_data_struct[2]
        _particles = domain_data_struct[4]

        self._map_to_grid(_gridx, _gridy, _particles, self._ibmf)

        return

    def force_flow(self, domain_data_struct):

        """
        Subroutine to compute immersed boundary forces
 
        Arguments
        ---------
        domain_data_struct : object list

        """

        _gridx = domain_data_struct[1]
        _gridy = domain_data_struct[2]
        _scalars = domain_data_struct[3]
        _particles = domain_data_struct[4]

        self._force_flow(_gridx, _gridy, _scalars, _particles, self._ibmf, self._velc)

        return
