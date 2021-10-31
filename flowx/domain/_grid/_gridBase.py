"""Module with implementation of the Grid classes."""

import numpy
import bubblebox.library.create as boxcreate

class GridBase(object):
    """Base class for the Grid."""

    type_ = 'base'

    def __init__(self, var_names, nx, ny, xmin, xmax, ymin, ymax,
                 user_bc_type=None, user_bc_val=None):
        """Initialize the Grid object and allocate the data.

        Parameters
        ----------
        var_names : list of strings
            List of names for the variables to create.
        nx : integer
            Number of cells in the x-direction.
        ny : integer
            Number of cells in the y-direction.
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
        # Domain information
        self.nx, self.ny = nx, ny  # number of cells in each direction
        self.xmin, self.xmax = xmin, xmax  # domain limits (x direction)
        self.ymin, self.ymax = ymin, ymax  # domain limits (y direction)
        self.dx = abs(self.xmax - self.xmin) / nx  # cell-width (x direction)
        self.dy = abs(self.ymax - self.ymin) / ny  # cell-width (y direction)
        self.set_gridline_coordinates()

        # Initialize the data
        self.num = len(var_names)
        self.vars = dict(zip(var_names, range(self.num)))
        self.data = boxcreate.Data()
        self.initialize_data()

        # Boundary condition information
        self.set_default_bc()
        if user_bc_type is not None and user_bc_val is not None:
            self.set_user_bc(user_bc_type, user_bc_val)
        self.fill_guard_cells(var_names)

    def __repr__(self):
        """Return a representation of the object."""
        return ("Grid:\n" +
                " - type: {}\n".format(type(self)) +
                " - size: {} x {}\n".format(self.nx, self.ny) +
                " - domain: [{}, {}] x [{}, {}]\n".format(self.xmin,
                                                          self.xmax,
                                                          self.ymin,
                                                          self.ymax) +
                " - number of variables: {}\n".format(self.num))

    def __del__(self):
        """Destructor"""
        self.data.purge()

    def set_gridline_coordinates(self):
        """Set the gridline coordinates.

        Method implemented in the child class.

        """
        raise NotImplementedError()

    def initialize_data(self):
        """Initialize the data with zeros.

        Method implemented in child classes.

        """
        raise NotImplementedError()

    def set_default_bc(self):
        """Set default boundary conditions (homogeneous Neumann)."""
        var_names = list(self.vars.keys())
        num = len(var_names)
        default_bc_type = 4 * ['neumann']
        self.bc_type = dict(zip(var_names, num * [default_bc_type]))
        default_bc_val = 4 * [0.0]
        self.bc_val = dict(zip(var_names, num * [default_bc_val]))

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

    def get_variable_indices(self, *var_names):
        """Get the grid index of given variable names.

        Parameters
        ----------
        var_names : tuple of strings
            The name of the variable(s).

        Returns
        -------
        indices : integer or list of integers
            Index of the grid variables.

        """
        indices = [self.vars[name] for name in var_names]
        return indices[0] if len(indices) == 1 else indices

    def set_values(self, var_name, values):
        """Set the values of a variable.

        Parameters
        ----------
        var_name : string
            Name of the variable.
        values : numpy.ndarray
            2D array with the values to set.

        """
        self.data[var_name][0,:,:,0] = values

    def get_values(self, var_name):
        """Get the data of a variable (as a copy).

        Parameters
        ----------
        var_name : string
            Name of the variable.

        Returns
        -------
        data : numpy.ndarray
            Variable's data as a 2D array of floats.

        """
        data = self.data[var_name][0,:,:,0]
        return data

    def set_value(self, var_name, i, j, value):
        """Set a value using (i, j) indexation.

        Parameters
        ----------
        var_name : string
            Name of the variable.
        i : integer
            Index in the x-direction.
        j : integer
            Index in the y-direction.
        value : float
            Value to set.

        """
        self.data[var_name][0,i,j,0] = value

    def get_value(self, var_name, i, j):
        """Get a value using (i, j) indexation.

        Parameters
        ----------
        var_name : string
            Name of the variable.
        i : integer
            Index in the x-direction.
        j : integer
            Index in the y-direction.

        Returns
        -------
        value : float
            Value of the variable at index (i, j).

        """
        value = self.data[var_name][0,i,j,0]
        return value

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
        self.data[eror][0,:,:,0] = numpy.abs(self.data[ivar][0,:,:,0] -
                                             self.data[asol][0,:,:,0])

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
        l2_norm = (numpy.sqrt(numpy.sum(self.data[eror][0,:,:,0]**2)) /
                   ((self.nx + 2) * (self.ny + 2)))

        return l2_norm

    def fill_guard_cells(self, var_names):
        """Fill value at guard cells for given variable names.

        Parameters
        ----------
        var_names : string or list of strings
            Name of variables to update.

        """
        # Convert single string to a list
        if type(var_names) is str:
            var_names = [var_names]
        # Set locations and delta value (either dx or dy)
        locs = ['left', 'right', 'bottom', 'top']
        deltas = [self.dx, self.dx, self.dy, self.dy]
        # Fill guard cells for each variable
        for name in var_names:
            bc_types = self.bc_type[name]
            bc_vals = self.bc_val[name]
            # Fill guard cells for each boundary
            for loc, delta, bc_type, bc_val in zip(locs, deltas,
                                                   bc_types, bc_vals):
                if bc_type == 'neumann':
                    self.fill_guard_cells_neumann(name, loc, bc_val, delta)
                elif bc_type == 'dirichlet':
                    self.fill_guard_cells_dirichlet(name, loc, bc_val)
                elif bc_type == 'outflow':
                    self.fill_guard_cells_dirichlet(name, loc, bc_val)
                elif bc_type == 'periodic':
                    self.fill_guard_cells_periodic(name, loc)
                elif bc_type == 'projection':
                    self.fill_guard_cells_projection(name, loc)
                elif bc_type == None:
                    None
                else:
                    raise ValueError('Boundary type "{}" not implemented'
                                     .format(bc_type))

    def fill_guard_cells_dirichlet(self, var_name, loc, bc_val):
        """Fill guard cells using a Dirichlet condition.

        Method implemented in child classes.

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
        raise NotImplementedError()

    def fill_guard_cells_neumann(self, var_name, loc, bc_val, delta):
        """Fill guard cells using a Neumann condition.

        Parameters
        ----------
        var_name : string
            Name of the variable to update.
        loc : string
            Boundary location;
            choices: ['left', 'right', 'bottom', 'top'].
        bc_val : float
            Neumann boundary value.
        delta : float
            Grid-cell width.

        """ 
        var = self.get_values(var_name)
        if loc == 'left':
            var[0, :] = bc_val * delta + var[1, :]
        elif loc == 'right':
            var[-1, :] = bc_val * delta + var[-2, :]
        elif loc == 'bottom':
            var[:, 0] = bc_val * delta + var[:, 1]
        elif loc == 'top':
            var[:, -1] = bc_val * delta + var[:, -2]
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))

    def fill_guard_cells_periodic(self, var_name, loc):
        """Fill guard cells with periodic BC.

        Parameters
        ----------
        var_name : string
            Name of the variable to update.
        loc : string
            Boundary location;
            choices: ['left', 'right', 'bottom', 'top'].
        """
        var = self.get_values(var_name)
        if loc == 'left':
            var[0, :] = var[-2, :]
        elif loc == 'right':
            var[-1, :] = var[1, :]
        elif loc == 'bottom':
            var[:, 0] = var[:, -2]
        elif loc == 'top':
            var[:, -1] = var[:, 1]
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))

    def fill_guard_cells_projection(self, var_name, loc):
        """Fill guard cells with projection BC.

        Parameters
        ----------
        var_name : string
            Name of the variable to update.
        loc : string
            Boundary location;
            choices: ['left', 'right', 'bottom', 'top'].
        """
        var = self.get_values(var_name)
        if loc == 'left':
            var[0, :] = 2*var[1, :] - var[2, :]
        elif loc == 'right':
            var[-1, :] = 2*var[-2, :] - var[-3, :]
        elif loc == 'bottom':
            var[:, 0] = 2*var[:, 1] - var[:, 2]
        elif loc == 'top':
            var[:, -1] = 2*var[:, -2] - var[:, -3]
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))
