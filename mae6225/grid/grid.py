"""Implementation of the class `Grid`."""

import numpy


class Grid(object):
    """A Grid object stores data about the grid and the variables."""

    def __init__(self, center_vars, nx, ny, xmin, xmax, ymin, ymax, user_bc_type, user_bc_val):
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

        user_bc_type: dictionary for user defined BCs

        user_bc_val: dictionary for user defined BC values

        """

        # Domain information
        self.nx, self.ny = nx, ny
        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.dx = abs(self.xmax - self.xmin) / nx
        self.dy = abs(self.ymax - self.ymin) / ny

        # Cell center co-ordinates
        self.x_center, self.y_center = self.get_cell_centered_coordinates()

        # Cell center data
        self.center_vars = center_vars
        self.num = len(center_vars)
        self.data = numpy.zeros((nx + 2, ny + 2, self.num), dtype=numpy.float64)

        # Boundary condition information
        self.bc_mask = numpy.full((self.num,1), False)
        self.bc_type,self.bc_val,self.bc_data_struct,self.bc_vars = self.set_bc()
        self.set_user_bc(user_bc_type,user_bc_val) 

    def __repr__(self):
        """Return a representation of the object."""
        return ("Grid:\n" +
                " - size: {} x {}\n".format(self.nx, self.ny) +
                " - domain: [{}, {}] x [{}, {}]\n".format(self.xmin,
                                                          self.xmax,
                                                          self.ymin,
                                                          self.ymax) +
                " - number of variables: {}\n".format(self.num))

    def set_bc(self):
        """Set boundary conditions for variables"""

        bc_type        = dict()
        bc_val         = dict()
        bc_vars        = dict()
        bc_data_struct = dict()

        var_count  = 0

        for key in self.center_vars:
            bc_type[key]        = ['neumann','neumann','neumann','neumann'] 
            bc_val[key]         = [0.,0.,0.,0.]
            bc_data_struct[key] = 'center' 
            bc_vars[key]        = var_count

            var_count = var_count + 1

        return bc_type,bc_val,bc_data_struct,bc_vars

    def set_user_bc(self,user_bc_type,user_bc_val):
        """Set user defined BC"""
 
        for key in user_bc_type:
            self.bc_type[key] = user_bc_type[key]
            self.bc_val[key]  = user_bc_val[key]

        return

    def set_values(self, idx, values):
        """Set the values of a variable.

        Parameters
        ----------
        idx : integer
            Index of the variable in the data array.
        values : numpy.ndarray
            2D array with the values to set.

        """
        self.data[:, :, self.center_vars[idx]] = values

    def get_values(self, idx):
        """Get the data of a variable (as a copy).

        Parameters
        ----------
        idx : integer
            Index of the variable in the data array.

        Returns
        -------
        data : numpy.ndarray
            Variable's data as a 2D array of floats.

        """
        data = self.data[:, :, self.center_vars[idx]]
        return data

    def set_value(self, idx, i, j, value):
        """Set a value using (i, j) indexation.

        Parameters
        ----------
        idx : integer
            Index of the variable in the data array.
        i : integer
            Index in the x-direction.
        j : integer
            Index in the y-direction.
        value : float
            Value to set.

        """
        self.data[i, j, self.center_vars[idx]] = value

    def get_value(self, idx, i, j):
        """Get a value using (i, j) indexation.

        Parameters
        ----------
        idx : integer
            Index of the variable in the data array.
        i : integer
            Index in the x-direction.
        j : integer
            Index in the y-direction.

        Returns
        -------
        value : float
            Value of the variable at index (i, j).

        """
        value = self.data[i, j, self.center_vars[idx]]
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

    def fill_guard_cells(self):

        """
        Function to apply boundary condition

        """
 
        for key in self.bc_vars:

            if(self.bc_mask[self.bc_vars[key]]):     
            
                if(self.bc_data_struct[key] == 'center'):

                    # xmin BC
                    if(self.bc_type[key][0] == 'neumann'):
                        self.data[:,0,self.center_vars[key]]  =  self.bc_val[key][0]*self.dx + self.data[:,1,self.center_vars[key]]
                    else:
                        self.data[:,0,self.center_vars[key]]  =  2.0*self.bc_val[key][0]     - self.data[:,1,self.center_vars[key]]
                    # xmax BC
                    if(self.bc_type[key][1] == 'neumann'):
                        self.data[:,-1,self.center_vars[key]] = -self.bc_val[key][1]*self.dx + self.data[:,-2,self.center_vars[key]]
                    else:
                        self.data[:,-1,self.center_vars[key]] =  2.0*self.bc_val[key][1]     - self.data[:,-2,self.center_vars[key]]

                    # ymin BC
                    if(self.bc_type[key][2] == 'neumann'):
                        self.data[0,:,self.center_vars[key]]  =  self.bc_val[key][2]*self.dy + self.data[1,:,self.center_vars[key]]
                    else:
                        self.data[0,:,self.center_vars[key]]  =  2.0*self.bc_val[key][2]     - self.data[1,:,self.center_vars[key]]
                    # ymax BC
                    if(self.bc_type[key][3] == 'neumann'):
                        self.data[-1,:,self.center_vars[key]] = -self.bc_val[key][3]*self.dy + self.data[-2,:,self.center_vars[key]]
                    else:
                        self.data[-1,:,self.center_vars[key]] = 2.0*self.bc_val[key][3]      - self.data[-2,:,self.center_vars[key]]
    
        return
