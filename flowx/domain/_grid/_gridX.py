"""Module with implementation of the Grid classes."""

import numpy

from . import GridBase

class GridFaceX(GridBase):
    """Class for a x-face grid."""

    type_ = 'x-face'

    def __init__(self,*args,**kwargs):
        """Call the constructor of the base class."""
        super().__init__(*args,**kwargs)

    @classmethod
    def check_gridtype(cls, gridtype):
        """Check if grid type if 'cell-centered'."""
        return gridtype == 'x-face'

    @staticmethod
    def initialize_data_attributes(nblocks,nxb,nyb,varlist):
        """Private method for initialization"""

        data_attributes = {'nblocks'   : nblocks,
                           'nxb'       : nxb+1,
                           'nyb'       : nyb,
                           'xguard'    : 0,
                           'yguard'    : 1,
                           'variables' : dict(zip(varlist,[None]*len(varlist)))}

        return data_attributes

    def set_gridline_coordinates(self):
        """Set the gridline coordinates."""
        for block in self.blocklist:
            block.x = numpy.linspace(block.xmin,
                                     block.xmax,
                                     num=self.nxb)
            block.y = numpy.linspace(block.ymin - block.dy / 2,
                                     block.ymax + block.dy / 2,
                                     num=self.nyb+2*self.yguard)
 
        self.x = numpy.linspace(self.xmin,
                                self.xmax,
                                num=self.nx+1)
        self.y = numpy.linspace(self.ymin - self.dy / 2,
                                self.ymax + self.dy / 2,
                                num=self.ny+2)

    @staticmethod
    def fill_guard_cells_dirichlet(blockdata, loc, bc_val):
        """Fill guard cells using a Dirichlet condition.

        Parameters
        ----------
        loc : string
            Boundary location;
            choices: ['left', 'right', 'bottom', 'top'].
        bc_val : float
            Neumann boundary value.

        """
        if loc == 'xlow':
            blockdata[:,:,0] = bc_val
        elif loc == 'xhigh':
            blockdata[:,:,-1] = bc_val
        elif loc == 'ylow':
            blockdata[:,0,:] = 2 * bc_val - blockdata[:,1,:]
        elif loc == 'yhigh':
            blockdata[:,-1,:] = 2 * bc_val - blockdata[:,-2,:]
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))
