"""Module with implementation of the Grid classes."""

import numpy

from . import GridBase

class GridFaceY(GridBase):
    """Class for a y-face grid."""

    type_ = 'y-face'

    def __init__(self,*args,**kwargs):
        """Call the constructor of the base class."""
        super().__init__(*args,**kwargs)

    @classmethod
    def check_gridtype(cls, gridtype):
        """Check if grid type if 'cell-centered'."""
        return gridtype == 'y-face'

    @staticmethod
    def initialize_data_attributes(nblocks,nxb,nyb,varlist):
        """Private method for initialization"""
        data_attributes = {'nblocks'   : nblocks,
                           'nxb'       : nxb+2,
                           'nyb'       : nyb+1,
                           'variables' : dict(zip(varlist,[None]*len(varlist)))}

        return data_attributes
 
    def set_gridline_coordinates(self):
        """Set the gridline coordinates."""
        for block in self.blocklist:
            block.x = numpy.linspace(block.xmin - block.dx / 2,
                                     block.xmax + block.dx / 2,
                                     num=self.nxb)
            block.y = numpy.linspace(block.ymin,
                                     block.ymax,
                                     num=self.nyb)
 
        self.x = numpy.linspace(self.xmin - self.dx / 2,
                                self.xmax + self.dx / 2,
                                num=self.nx+2)
        self.y = numpy.linspace(self.ymin,
                                self.ymax,
                                num=self.ny+1)
  
    @staticmethod
    def fill_guard_cells_dirichlet(blockdata, loc, bc_val):
        """Fill guard cells using a Dirichlet condition.

        Parameters
        ----------
            Name of the variable to update.
        loc : string
            Boundary location;
            choices: ['left', 'right', 'bottom', 'top'].
        bc_val : float
            Neumann boundary value.

        """
        if loc == 'xlow':
            blockdata[:,:,0] = 2 * bc_val - blockdata[:,:,1]
        elif loc == 'xhigh':
            blockdata[:,:,-1] = 2 * bc_val - blockdata[:,:,-2]
        elif loc == 'ylow':
            blockdata[:,0,:] = bc_val
        elif loc == 'yhigh':
            blockdata[:,-1,:] = bc_val
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))
