"""Module with implementation of the Grid classes."""

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

    def fill_guard_cells_dirichlet(block, ivar, loc, bc_val):
        """Fill guard cells using a Dirichlet condition.

        Parameters
        ----------
        block : BubbleBox Block object
        ivar : string
            Name of the variable to update.
        loc : string
            Boundary location;
            choices: ['left', 'right', 'bottom', 'top'].
        bc_val : float
            Neumann boundary value.

        """
        if loc == 'left':
            block[ivar][:,:,0]  = 2 * bc_val - block[ivar][:,:,1]
        elif loc == 'right':
            block[ivar][:,:,-1] = 2 * bc_val - block[ivar][:,:,-2]
        elif loc == 'bottom':
            block[ivar][:,0,:]  = 2 * bc_val - block[ivar][:,1,:]
        elif loc == 'top':
            block[ivar][:,-1,:] = 2 * bc_val - block[ivar][:,-2,:]
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))
