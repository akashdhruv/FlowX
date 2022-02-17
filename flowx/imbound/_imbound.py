"""Implementation of the immersed boundary module"""

import time

from . import _interface


class ImBound(object):
    def __init__(
        self, domain_data_struct=[None] * 5, imbound_vars=[None] * 4, imbound_info=None
    ):
        """
        Constructor for the imbound unit

        Arguments
        ---------
        domain_data_struct : object list
                             [gridc, gridx, gridy, scalars, particles]

        imbound_vars       : list
                             List of string for field variables required by ins unit
                             imbound_vars[0] --> indicator variable for immersed boundary
                             imbound_vars[1] --> velocity variable
                             imbound_vars[2] --> dynamic x grid variable
                             imbound_vars[3] --> dynamic y grid variable

        imbound_info : Dictionary of keyword arguments

        'with_ib' - keyword to indicate if immersed boundary is present or not
                  - True if present
                  - False if not present

        imbound_info['with_ib'] = False --> default
                                = True

        Error generated if kwargs['with_ib'] is True and imbound_vars is empty

        """
        # ----------Create images for other unit objects and variables----------------------------
        (
            self._gridc,
            self._gridx,
            self._gridy,
            self._scalars,
            self._particles,
        ) = domain_data_struct

        imbound_vars.extend([None] * (4 - len(imbound_vars)))

        self._ibmf, self._velc, self._ibmx, self._ibmy = imbound_vars

        # --------------Set default parameters for the current unit--------------------------------
        self._options = {
            "with_ib": False,
            "mapping_type": "ann",
            "verbose": False,
            "ntrees": 20,
            "nquery_trees": 2,
            "nquery_trace": 30,
            "ib_type": "rigid",
            "lset_redistance": 3,
            "extrap_solid": 20,
            "monitor": False,
            "nthreads": 1,
            "backend": "serial",
        }

        self._mapping_type = {
            "classical": _interface.utils.classical_search,
            "ann": _interface.utils.ann_search,
            "shapely": _interface.utils.shapely_search,
        }

        self._mapping_time = 0.0
        self._mapping_ites = 0
        self._advection_time = 0.0

        # ---------------------------Read user parameters---------------------------------------
        if imbound_info:
            for key in imbound_info:
                self._options[key] = imbound_info[key]

        # --------------------------Setup unit based on user/default parameters---------------
        if (
            not self._options["with_ib"]
            or all(var is None for var in imbound_vars)
            or None in domain_data_struct
        ):
            self._force_flow = _interface.stub.force_flow
            self._map_to_grid = _interface.stub.map_to_grid
            self._advect = _interface.stub.advect
            print(
                "Warning: Immersed Boundary unit is a stub, no forcing will be applied."
            )

        else:
            if self._options["ib_type"] == "rigid":
                self._force_flow = _interface.rigid.force_flow
                self._map_to_grid = _interface.rigid.map_to_grid
                self._advect = _interface.stub.advect

            elif self._options["ib_type"] == "visco":
                self._force_flow = _interface.visco.force_flow
                self._map_to_grid = _interface.visco.map_to_grid
                self._advect = _interface.visco.advect

            self._search_function = self._mapping_type[self._options["mapping_type"]]

        return

    def map_to_grid(self):
        """
        Subroutine to map immersed boundary on grid
        """
        t1 = time.time()
        self._mapping_ites = self._map_to_grid(
            self._gridc,
            self._particles,
            self._ibmf,
            self._ibmx,
            self._ibmy,
            self._search_function,
            self._options,
        )
        t2 = time.time()
        self._mapping_time = t2 - t1

    def force_flow(self):
        """
        Subroutine to compute immersed boundary forces
        """
        self._force_flow(
            self._gridc,
            self._gridx,
            self._gridy,
            self._scalars,
            self._particles,
            self._ibmf,
            self._ibmx,
            self._ibmy,
            self._velc,
            self._options,
        )

    def advect(self):
        """
        Subroutine to compute immersed boundary forces
        """
        t1 = time.time()
        self._advect(
            self._gridc,
            self._gridx,
            self._gridy,
            self._scalars,
            self._ibmf,
            self._ibmx,
            self._ibmy,
            self._velc,
            self._options,
        )
        t2 = time.time()

        self._advection_time = t2 - t1
