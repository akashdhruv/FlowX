"""Module with implementation of the Grid classes."""

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

    def _initialize_data_attributes(self,xblocks,yblocks,nxb,nyb,varlist):
        """Private method for initialization"""

        data_attributes = {'nblocks'   : xblocks*yblocks,
                           'nxb'       : nxb+1,
                           'nyb'       : nyb+2,
                           'variables' : dict(zip(varlist,[None]*len(varlist)))}

        return data_attributes
