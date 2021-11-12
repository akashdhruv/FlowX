"""Module with implementation of the Grid classes."""

from bubblebox.library.create import Dataset,Block,Data
from bubblebox.library.utilities import Action
import numpy
import pymorton

class GridBase(Dataset):
    """Base class for the Grid."""

    type_ = 'base'

    def __init__(self, varlist, nx, ny, xmin, xmax, ymin, ymax, xblocks=1, yblocks=1,
                                                    user_bc_type=None, user_bc_val=None):
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

        # Organize data at coarsest level
        dx,dy = abs(xmax-xmin)/nx,abs(ymax-ymin)/ny
        nxb,nyb = int(nx/xblocks),int(ny/yblocks)
        levels = None

        # Save grid attributes at coarsest level
        self.nx,self.ny = nx,ny
        self.dx,self.dy = dx,dy

        # Initialize block attributes
        block_attributes = self.__class__.initialize_block_attributes(xblocks,yblocks,dx,dy,
                                                                      xmin,xmax,ymin,ymax,levels)

        # Initialize data attributes
        nblocks = len(block_attributes)
        data_attributes = self.__class__.initialize_data_attributes(nblocks,nxb,nyb,varlist)

        # Create data and block objects
        data = Data(**data_attributes)
        blocklist = [Block(data,**attributes) for attributes in block_attributes]

        # Call base class constructor
        super().__init__(blocklist,data)
       
        # Set gridline coordinates
        self.set_gridline_coordinates()

        # Set boundary blocks
        #self.set_domain_boundaries()

        # Boundary condition information
        self.bc_type = {}
        self.bc_val = {}

        self.set_default_bc(varlist)
        if user_bc_type is not None and user_bc_val is not None:
            self.set_user_bc(user_bc_type, user_bc_val)

        self.fill_guard_cells(varlist)

    def __del__(self):
        """Destructor"""
        self.purge()

    @staticmethod
    def initialize_block_attributes(xblocks,yblocks,dx,dy,xmin,xmax,ymin,ymax,levels):
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

    @staticmethod
    def initialize_data_attributes(nblocks,nxb,nyb,varlist):
        """Private method for initialization"""
        raise NotImplementedError

    def set_gridline_coordinates(self):
        """Set the gridline coordinates."""
        raise NotImplementedError

    def addvar(self,varkey):
        """Add a variable"""
        super().addvar(varkey)
        self.set_default_bc(varkey)

    def set_default_bc(self,varlist):
        """Set default boundary conditions (homogeneous Neumann)."""
        if type(varlist) is str:
            varlist = [varlist]

        default_bc_type = 4 * ['neumann']
        default_bc_val = 4 * [0.0]   

        num = len(varlist)
 
        self.bc_type = {**self.bc_type, **dict(zip(varlist, num * [default_bc_type]))}
        self.bc_val = {**self.bc_val, **dict(zip(varlist, num * [default_bc_val]))}

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

    def compute_error(self, eror, ivar, asol):
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
            block[eror] = numpy.abs(block[ivar] - block[asol])

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
        l2_norm = 0.

        for block in self.blocklist:
            l2_norm = l2_norm + (numpy.sqrt(numpy.sum(block[eror]**2)) /
                      ((self.nxb+2*self.xguard) * (self.nyb+2*self.yguard)))

        return l2_norm/self.nblocks
 
    def fill_guard_cells(self, varlist, **kwargs):
        """Fill value at guard cells for given variable names.

        Parameters
        ----------
        varlist : string or list of strings
            Name of variables to update.

        """
        self.halo_exchange(varlist, **kwargs)

        # Convert single string to a list
        if type(varlist) is str:
            varlist = [varlist]

        locations = ['xlow','xhigh','ylow','yhigh']

        # TODO add a call to exchange data between blocks
        # TODO figure out how to tag blocks at boundary etc.
        # TODO make this efficient
        for varkey in varlist:
            bc_type_var = self.bc_type[varkey]
            bc_val_var = self.bc_val[varkey]
               
            for block in self.blocklist:
                deltas = [block.dx, block.dx, block.dy, block.dy]
                neighbors = [block.neighdict[location] for location in locations]
                blockdata = block[varkey]

                for location,neighbor,delta,bc_type,bc_val in zip(locations,neighbors,deltas,
                                                                  bc_type_var,bc_val_var):
                    if neighbor is None:
                        if bc_type == 'neumann':
                            self.__class__.fill_guard_cells_neumann(blockdata,location,bc_val,delta)
                        elif bc_type == 'dirichlet':
                            self.__class__.fill_guard_cells_dirichlet(blockdata,location,bc_val)
                        elif bc_type == 'outflow':
                            self.__class__.fill_guard_cells_dirichlet(blockdata,location,bc_val)
                        elif bc_type == 'projection':
                            self.__class__.fill_guard_cells_projection(blockdata,location)
                        elif bc_type == None:
                            None
                        else:
                            raise ValueError('Boundary type "{}" not implemented'.format(bc_type))

    @staticmethod
    def fill_guard_cells_dirichlet(blockdata, loc, bc_val):
        """Fill guard cells using a Dirichlet condition.

        Method implemented in child classes.

        Parameters
        ----------
        loc : string
            Boundary location;
            choices: ['left', 'right', 'bottom', 'top'].
        bc_val : float
            Neumann boundary value.

        """
        raise NotImplementedError()

    @staticmethod
    def fill_guard_cells_neumann(blockdata, loc, bc_val, delta):
        """Fill guard cells using a Neumann condition.

        Parameters
        ----------
        loc : string
            Boundary location;
            choices: ['left', 'right', 'bottom', 'top'].
        bc_val : float
            Neumann boundary value.
        delta : float
            Grid-cell width.

        """ 
        if loc == 'xlow':
            blockdata[:,:,0]  = bc_val * delta + blockdata[:,:,1]
        elif loc == 'xhigh':
            blockdata[:,:,-1] = bc_val * delta + blockdata[:,:,-2]
        elif loc == 'ylow':
            blockdata[:,0,:]  = bc_val * delta + blockdata[:,1,:]
        elif loc == 'yhigh':
            blockdata[:,-1,:] = bc_val * delta + blockdata[:,-2,:]
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))

    @staticmethod
    def fill_guard_cells_projection(blockdata, loc):
        """Fill guard cells with projection BC.

        Parameters
        ----------
        loc : string
            Boundary location;
            choices: ['left', 'right', 'bottom', 'top'].
        """
        if loc == 'xlow':
            blockdata[:,:,0]  = 2*blockdata[:,:,1]  - blockdata[:,:,2]
        elif loc == 'xhigh':
            blockdata[:,:,-1] = 2*blockdata[:,:,-2] - blockdata[:,:,-3]
        elif loc == 'ylow':
            blockdata[:,0,:]  = 2*blockdata[:,1,:]  - blockdata[:,2,:]
        elif loc == 'yhigh':
            blockdata[:,-1,:] = 2*blockdata[:,-2,:] - blockdata[:,-3,:]
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))
