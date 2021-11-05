"""Module with implementation of the Grid classes."""

import numpy

from . import GridBase

class GridCellCentered(GridBase):
    """Class for a cell-centered grid."""

    type_ = 'cell-centered'

    def __init__(self,*args,**kwargs):
        """Call the constructor of the base class."""
        super().__init__(*args,**kwargs)

    @classmethod
    def check_gridtype(cls, gridtype):
        """Check if grid type if 'cell-centered'."""
        return gridtype == 'cell-centered'

    def _initialize_data_attributes(self,xblocks,yblocks,nxb,nyb,varlist):
        """Private method for initialization"""

        data_attributes = {'nblocks'   : xblocks*yblocks,
                           'nxb'       : nxb+2,
                           'nyb'       : nyb+2,
                           'variables' : dict(zip(varlist,[None]*len(varlist)))}

        return data_attributes
 
    def set_gridline_coordinates(self):
        """Set the gridline coordinates."""
        for block in self.blocklist:
            block.x = numpy.linspace(block.xmin - block.dx / 2,
                                     block.xmax + block.dx / 2,
                                     num=self.nxb)
            block.y = numpy.linspace(block.ymin - block.dy / 2,
                                     block.ymax + block.dy / 2,
                                     num=self.nyb)
  
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
            blockdata[:,:,0]  = 2 * bc_val - blockdata[:,:,1]
        elif loc == 'xhigh':
            blockdata[:,:,-1] = 2 * bc_val - blockdata[:,:,-2]
        elif loc == 'ylow':
            blockdata[:,0,:]  = 2 * bc_val - blockdata[:,1,:]
        elif loc == 'yhigh':
            blockdata[:,-1,:] = 2 * bc_val - blockdata[:,-2,:]
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))
