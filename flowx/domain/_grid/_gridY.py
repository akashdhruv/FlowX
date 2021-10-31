"""Module with implementation of the Grid classes."""

import numpy
from . import GridBase

class GridFaceY(GridBase):
    """Class for a y-face centered grid."""

    type_ = 'y-face'

    def __init__(self, *args, **kwargs):
        """Call the constructor of the base class."""
        super(GridFaceY, self).__init__(*args, **kwargs)

    @classmethod
    def check_gridtype(cls, gridtype):
        """Check if grid type if 'y-face'."""
        return gridtype == 'y-face'

    def set_gridline_coordinates(self):
        """Set the gridline coordinates."""
        self.x = numpy.linspace(self.xmin - self.dx / 2,
                                self.xmax + self.dx / 2,
                                num=self.nx + 2)
        self.y = numpy.linspace(self.ymin, self.ymax, num=self.ny + 1)

    def initialize_data(self):
        """Initialize the data with zeros."""
        self.data.nxb = self.nx+2
        self.data.nyb = self.ny+1
        [self.data.addvar(var) for var in self.vars]

    def fill_guard_cells_dirichlet(self, var_name, loc, bc_val):
        """Fill guard cells using a Dirichlet condition.

        Parameters
        ----------
        var_name : string
            Name of the variable to update.
        loc : string
            Boundary location;
            choices: ['left', 'right', 'bottom', 'top'].
        bc_val : float
            Neumann boundary value.

        """
        var = self.get_values(var_name)
        if loc == 'left':
            var[0, :] = 2 * bc_val - var[1, :]
        elif loc == 'right':
            var[-1, :] = 2 * bc_val - var[-2, :]
        elif loc == 'bottom':
            var[:, 0] = bc_val
        elif loc == 'top':
            var[:, -1] = bc_val
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))
