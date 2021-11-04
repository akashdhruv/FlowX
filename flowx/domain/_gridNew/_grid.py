"""Module with implementation of the Grid classes."""

from . import GridBase
from . import GridFaceX
from . import GridFaceY
from . import GridCellCentered

def Grid(gridtype, *args, **kwargs):
    """Return an instance of the GridBase child class based on type.

    Parameters
    ----------
    gridtype : string
        Type of grid;
        choices: ['cell-centered', 'x-face', 'y-face'].
    varlist : list of strings
            List of names for the variables to create.
    nx : integer
        Number of cells in the x-direction.
    ny : integer
        Number of cells in the y-direction.
    xblocks : integer
        Number of blocks in the x-direction
    yblocks : integer
        Number of blocks in the y-direction
    xmin : float
        Domain limit at the left side.
    xmax : float
        Domain limit at the right side.
    ymin : float
        Domain limit at the bottom side.
    ymax : float
        Domain limit at the top side.
    user_bc_type : dictionary of (string, list) items
        User-defined boundary types to overwrite default ones.
    user_bc_val : dictionary of (string, list) items
        User-defined boundary values to overwrite default ones.

    Returns
    -------
    obj : instance of the GridBase child
        The grid object.

    """
    for cls in GridBase.__subclasses__():
        if cls.check_gridtype(gridtype):
            return cls(*args, **kwargs)
    raise ValueError('Parameter "gridtype" should be either '
                     '"cell-centered", "x-face", or "y-face"')
