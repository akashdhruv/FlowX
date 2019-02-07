"""Module with the implementation of the class `Grid`."""

import numpy


class Grid(object):
    """A Grid object stores data about the grid and the variables."""

    def __init__(self, center_vars, nx, ny, xmin, xmax, ymin, ymax,
                 user_bc_type=None, user_bc_val=None):
        """Initialize the Grid object and allocate the data.

        Parameters
        ----------
        center_vars : dictionary for cell center variables

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
        self.nx, self.ny = nx, ny
        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.dx = abs(self.xmax - self.xmin) / nx
        self.dy = abs(self.ymax - self.ymin) / ny

        # Cell-centered coordinates
        self.x_center, self.y_center = self.get_cell_centered_coordinates()

        # Cell-centered data
        self.num = len(center_vars)
        self.center_vars = dict(zip(center_vars, range(self.num)))
        self.data = numpy.zeros((nx + 2, ny + 2, self.num))

        # Boundary condition information
        self.set_default_bc()
        if user_bc_type is not None and user_bc_val is not None:
            self.set_user_bc(user_bc_type, user_bc_val)

    def __repr__(self):
        """Return a representation of the object."""
        return ("Grid:\n" +
                " - size: {} x {}\n".format(self.nx, self.ny) +
                " - domain: [{}, {}] x [{}, {}]\n".format(self.xmin,
                                                          self.xmax,
                                                          self.ymin,
                                                          self.ymax) +
                " - number of variables: {}\n".format(self.num))

    def set_default_bc(self):
        """Set default boundary conditions (homogeneous Neumann)."""
        var_names = list(self.center_vars.keys())
        num = len(var_names)
        self.bc_type = dict(zip(var_names, num * 4 * ['neumann']))
        self.bc_val = dict(zip(var_names, num * 4 * [0.0]))
        self.bc_data_struct = dict(zip(var_names, num * ['center']))

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
        # Oeverwrite default boundary values
        self.bc_val = {**self.bc_val, **user_bc_val}

    def get_variable_indices(self, var_names):
        """Get the grid index of given variable names.

        Parameters
        ----------
        var_names : string or list of strings
            The name of the variable(s).

        Returns
        -------
        indices : list of integers
            Index of the grid variables.

        """
        # Convert single string to list
        if var_names is str:
            var_names = [var_names]
        indices = []
        for name in var_names:
            indices.append(self.center_vars[name])
        return indices

    def set_values(self, var_name, values):
        """Set the values of a variable.

        Parameters
        ----------
        var_name : string
            Name of the variable.
        values : numpy.ndarray
            2D array with the values to set.

        """
        idx = self.center_vars[var_name]
        self.data[:, :, idx] = values

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
        idx = self.center_vars[var_name]
        data = self.data[:, :, idx]
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
        idx = self.center_vars[var_name]
        self.data[i, j, idx] = value

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
        idx = self.center_vars[var_name]
        value = self.data[i, j, idx]
        return value

    def get_cell_centered_coordinates(self):
        """Return the cell-centered gridline coordinates.

        The gridline coordinates also contain the coordinate of the
        boundary ghost cell.

        Returns
        -------
        x : numpy.ndarray
            x-coordinates along a gridline as a 1D array of floats.
        y : numpy.ndarray
            y-coordinates along a gridline as a 1D array of floats.

        """
        x = numpy.linspace(self.xmin - self.dx / 2, self.xmax + self.dx / 2,
                           num=self.nx + 2)
        y = numpy.linspace(self.ymin - self.dy / 2, self.ymax + self.dy / 2,
                           num=self.ny + 2)
        return x, y

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
            # Only fill guard cells for cell-centered variables
            if self.bc_data_struct[name] == 'center':
                bc_types = self.bc_type[name]
                bc_vals = self.bc_val[name]
                # Fill guard cells for each boundary
                for loc, delta, bc_type, bc_val in zip(locs, deltas,
                                                       bc_types, bc_vals):
                    if bc_type == 'neumann':
                        self.fill_guard_cells_neumann(name, loc, bc_val, delta)
                    elif bc_type == 'dirichlet':
                        self.fill_guard_cells_dirichlet(name, loc, bc_val)
                    else:
                        raise ValueError('Boundary type "{}" not implemented'
                                         .format(bc_type))

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
        idx = self.center_vars[var_name]
        if loc == 'left':
            self.data[:, 0, idx] = bc_val * delta + self.data[:, 1, idx]
        elif loc == 'right':
            self.data[:, -1, idx] = bc_val * delta + self.data[:, -2, idx]
        elif loc == 'bottom':
            self.data[0, :, idx] = bc_val * delta + self.data[1, :, idx]
        elif loc == 'top':
            self.data[-1, :, idx] = bc_val * delta + self.data[-2, :, idx]
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))

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
        idx = self.center_vars[var_name]
        if loc == 'left':
            self.data[:, 0, idx] = 2.0 * bc_val - self.data[:, 1, idx]
        elif loc == 'right':
            self.data[:, -1, idx] = 2.0 * bc_val - self.data[:, -2, idx]
        elif loc == 'bottom':
            self.data[0, :, idx] = 2.0 * bc_val - self.data[1, :, idx]
        elif loc == 'top':
            self.data[-1, :, idx] = 2.0 * bc_val - self.data[-2, :, idx]
        else:
            raise ValueError('Unknown boundary location "{}"'.format(loc))
