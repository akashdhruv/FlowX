"""Module with implementation of the Grid classes."""

import bubblebox.library as library
import pymorton

class GridBase(library.create.Dataset):
    """Base class for the Grid."""

    type_ = 'base'

    def __init__(self, varlist, nx, ny, xblocks, yblocks, 
                 xmin, xmax, ymin, ymax, user_bc_type=None, user_bc_val=None):
        """
        Initialize the Grid object and allocate the data.

        Parameters
        ----------
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
        """ 
        # Perform checks
        if nx%xblocks or ny%yblocks:
            raise ValueError('[flowx.domain.GridBase]:(nx,ny) must be exactly '+
                               'divisible by (xblocks,yblocks)')
        elif (xblocks%2 or yblocks%2) and xblocks!=1 and yblocks!=1:
            raise ValueError('[flowx.domain.GridBase]:(xblocks,yblocks) must be exactly '+
                               'divisible by 2')

        # Organize data
        nxb,nyb = [int(nx/xblocks),int(ny/yblocks)]
        dx,dy   = abs(xmax-xmin)/nx,abs(ymax-ymin)/ny

        # Initialize attributes
        block_attributes = self._initialize_block_attributes(xblocks,yblocks,dx,dy,xmin,xmax,ymin,ymax)
        data_attributes = self._initialize_data_attributes(xblocks,yblocks,nxb,nyb,varlist)

        # Create data and block objects
        data = library.create.Data(**data_attributes)
        blocklist = [library.create.Block(data,**attributes) for attributes in block_attributes]

        super().__init__(blocklist,data)

    def __del__(self):
        """Destructor"""
        self.purge()

    def _initialize_block_attributes(self,xblocks,yblocks,dx,dy,xmin,xmax,ymin,ymax):
        """Private method for initialization"""
        block_attributes = []

        for lblock in range(xblocks*yblocks):

            iloc,jloc = pymorton.deinterleave2(lblock)
            imin,imax = [xmin + (iloc/xblocks)*(xmax-xmin), xmin + ((iloc+1)/xblocks)*(xmax-xmin)]
            jmin,jmax = [ymin + (jloc/yblocks)*(ymax-ymin), ymin + ((jloc+1)/yblocks)*(ymax-ymin)]

            block_attributes.append({'dx'   : dx,
                                     'dy'   : dy,
                                     'xmin' : imin,
                                     'xmax' : imax,
                                     'ymin' : jmin,
                                     'ymax' : jmax,
                                     'tag'  : lblock})
        return block_attributes

    def _initialize_data_attributes(self,xblocks,yblocks,nxb,nyb,varlist):
        """Private method for initialization"""
        raise NotImplementedError

    def set_default_bc(self):
        """Set default boundary conditions (homogeneous Neumann)."""
        num = len(self.varlist)
        default_bc_type = 4 * ['neumann']
        self.bc_type = dict(zip(self.varlist, num * [default_bc_type]))
        default_bc_val = 4 * [0.0]
        self.bc_val = dict(zip(self.varlist, num * [default_bc_val]))

    def set_user_bc(self, user_bc_type, user_bc_val):
        """Overwrite default boundary conditions with user-provided ones.

        Parameters
        ----------
        user_bc_type : dictionary of (string, list) items
            User-defined boundary types.
        user_bc_val : dictionary of (string, list) items
            User-defined boundary values.

        """
        # Overwrite default boundary types
        self.bc_type = {**self.bc_type, **user_bc_type}
        # Overwrite default boundary values
        self.bc_val = {**self.bc_val, **user_bc_val}

    def update_bc_val(self, user_bc_val):
        """Overwrite boundary condition values with user-provided ones.

        Parameters
        ----------
        user_bc_val : dictionary of (string, list) items
            User-defined boundary values.

        """
        self.bc_val = {**self.bc_val, **user_bc_val}
 
    def update_bc_type(self, user_bc_type):

        self.bc_type = {**self.bc_type, **user_bc_type}

    def get_error(self, eror, ivar, asol):
        """Compute the error between the numerical and analytical solutions.

        Error is defined as the absolute difference between the two solutions.

        Arguments
        ---------
        eror : string
            Name of the grid variable of the error.
        ivar : string
            Name of the grid variable of the numerical solution.
        asol : string
            Name of the grid variable of the analytical solution.

        """
        for block in self.blocklist:
            block[eror][0,:,:] = numpy.abs(block[ivar][0,:,:] -
                                           block[asol][0,:,:])

    def get_l2_norm(self, eror):
        """Compute the L2 norm for a given variable.

        Arguments
        ---------
        eror : string
            Name of the grid variable for which norm is desired

        Returns
        -------
        l2_norm : float
            The L2-norm.

        """
        # TODO add treatment for multiple blocks
        for block in self.blocklist:
            l2_norm = (numpy.sqrt(numpy.sum(block[eror][0,:,:]**2)) /
                      ((self.nx + 2) * (self.ny + 2)))

        return l2_norm
 
    def fill_guard_cells(self, varlist):
        """Fill value at guard cells for given variable names.

        Parameters
        ----------
        varlist : string or list of strings
            Name of variables to update.

        """
        # Convert single string to a list
        if type(varlist) is str:
            varlist = [varlist]
        # Set locations and delta value (either dx or dy)
        locs = ['left', 'right', 'bottom', 'top']

        # TODO add a call to exchange data between blocks
        # TODO figure out how to tag blocks at boundary etc.
        for block in self.blocklist:
            deltas = [block.dx, block.dx, block.dy, block.dy]
            # Fill guard cells for each variable
            for ivar in varlist:
                bc_types = self.bc_type[ivar]
                bc_vals = self.bc_val[ivar]
                # Fill guard cells for each boundary
                for loc, delta, bc_type, bc_val in zip(locs, deltas,
                                                       bc_types, bc_vals):

                    if bc_type == 'neumann':
                        fill_guard_cells_neumann(block,ivar,loc,bc_val,delta)
                    elif bc_type == 'dirichlet':
                        fill_guard_cells_dirichlet(block,ivar,loc,bc_val)
                    elif bc_type == 'outflow':
                        fill_guard_cells_dirichlet(block,ivar,loc,bc_val)
                    elif bc_type == 'periodic':
                        fill_guard_cells_periodic(block,ivar,loc)
                    elif bc_type == 'projection':
                        fill_guard_cells_projection(block,ivar,loc)
                    elif bc_type == None:
                        None
                    else:
                        raise ValueError('Boundary type "{}" not implemented'
                                        .format(bc_type))

    @staticmethod
    def fill_guard_cells_dirichlet(block, ivar, loc, bc_val):
        """Fill guard cells using a Dirichlet condition.

        Method implemented in child classes.

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
        raise NotImplementedError()

    @staticmethod
    def fill_guard_cells_neumann(block, ivar, loc, bc_val, delta):
        """Fill guard cells using a Neumann condition.

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
        delta : float
            Grid-cell width.

        """ 
        if loc == 'left':
            block[ivar][:,:,0]  = bc_val * delta + block[ivar][:,:,1]
        elif loc == 'right':
            block[ivar][:,:,-1] = bc_val * delta + block[ivar][:,:,-2]
        elif loc == 'bottom':
            block[ivar][:,0,:]  = bc_val * delta + block[ivar][:,1,:]
        elif loc == 'top':
            block[ivar][:,-1,:] = bc_val * delta + block[ivar][:,-2,:]
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))

    @staticmethod
    def fill_guard_cells_periodic(block, ivar, loc):
        """Fill guard cells with periodic BC.

        Parameters
        ----------
        block : BubbleBox Block object
        ivar : string
            Name of the variable to update.
        loc : string
            Boundary location;
            choices: ['left', 'right', 'bottom', 'top'].
        """
        if loc == 'left':
            block[ivar][:,:,0]  = block[ivar][:,:,-2]
        elif loc == 'right':
            block[ivar][:,:,-1] = block[ivar][:,:,1]
        elif loc == 'bottom':
            block[ivar][:,0,:]  = block[ivar][:,-2,:]
        elif loc == 'top':
            block[ivar][:,-1,:] = block[ivar][:,1,:]
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))

    @staticmethod
    def fill_guard_cells_projection(block, ivar, loc):
        """Fill guard cells with projection BC.

        Parameters
        ----------
        block : BubbleBox Block object
        ivar : string
            Name of the variable to update.
        loc : string
            Boundary location;
            choices: ['left', 'right', 'bottom', 'top'].
        """
        if loc == 'left':
            block[ivar][:,:,0]  = 2*block[ivar][:,:,1]  - block[ivar][:,:,2]
        elif loc == 'right':
            block[ivar][:,:,-1] = 2*block[ivar][:,:,-2] - block[ivar][:,:,-3]
        elif loc == 'bottom':
            block[ivar][:,0,:]  = 2*block[ivar][:,1,:]  - block[ivar][:,2,:]
        elif loc == 'top':
            block[ivar][:,-1,:] = 2*block[ivar][:,-2,:] - block[ivar][:,-3,:]
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))
