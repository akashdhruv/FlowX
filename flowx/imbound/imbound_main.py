"""Implementation of the immersed boundary module"""

from flowx.imbound.imbound_interface import imbound_interface

class imbound_main(imbound_interface):

    def __init__(self, domain_data_struct=[None]*5, imbound_vars=[None]*2, imbound_info=None):

        """
        Constructor for the imbound unit

        Arguments
        ---------

        domain_data_struct : object list
              [gridc, gridx, gridy, scalars, particles]

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
        from flowx.imbound.solvers.search_algorithms import classical_search, ann_search, shapely_search

        self._gridc, self._gridx, self._gridy, self._scalars, self._particles = domain_data_struct

        self._ibmf, self._velc = imbound_vars

        self._options = {'with_ib' : False, 'mapping_type' : 'ann', 'verbose' : False, \
                         'ntrees' : 20, 'nquery_trees' : 2, 'nquery_trace' : 30}

        self._force_flow  = force_flow_stub
        self._map_to_grid = map_to_grid_stub

        self._mapping_time = 0.0
        self._mapping_ites = 0

        if imbound_info:
            for key in imbound_info: self._options[key] = imbound_info[key]

        if self._options['with_ib']:
            self._force_flow = force_flow_levelset
            self._map_to_grid = map_to_grid_levelset

        self._mapping_type = {'classical': classical_search, 'ann': ann_search, 'shapely' : shapely_search}

        self._search_function = self._mapping_type[self._options['mapping_type']]

        if self._options['with_ib'] and None in imbound_vars: 
            raise ValueError('imbound_vars cannot be empty when body flag is true')

        if not self._options['with_ib']:
             print('Warning: Immersed Boundary unit is a stub, no forcing will be applied.')
 
        return
 
    def map_to_grid(self):
        """
        Subroutine to map immersed boundary on grid
 
        """
        import time

        t1 = time.time()
        self._mapping_ites = self._map_to_grid(self._gridc, self._gridx, self._gridy, \
                                               self._particles, self._ibmf, self._search_function, self._options)
        t2 = time.time()

        self._mapping_time = t2-t1

        return

    def force_flow(self):

        """
        Subroutine to compute immersed boundary forces

        """

        self._force_flow(self._gridx, self._gridy, self._scalars, self._particles, self._ibmf, self._velc)

        return
