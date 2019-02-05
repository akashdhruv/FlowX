"""Implementation of the class `Grid`."""

import numpy


class Grid(object):
    """A Grid object stores data about the grid and the variables."""

    def __init__(self, num, nx, ny, xmin, xmax, ymin, ymax):
        """Initialize the Grid object and allocate the data.

        Parameters
        ----------
        num : integer
            Number of variables on the grid.
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

        """
        self.num = num
        self.nx, self.ny = nx, ny
        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.dx = abs(self.xmax - self.xmin) / nx
        self.dy = abs(self.ymax - self.ymin) / ny
        self.x_center, self.y_center = self.get_cell_centered_coordinates()
        self.data = numpy.zeros((nx + 2, ny + 2, num), dtype=numpy.float64)

    def __repr__(self):
        """Return a representation of the object."""
        return ("Grid:\n" +
                " - size: {} x {}\n".format(self.nx, self.ny) +
                " - domain: [{}, {}] x [{}, {}]\n".format(self.xmin,
                                                          self.xmax,
                                                          self.ymin,
                                                          self.ymax) +
                " - number of variables: {}\n".format(self.num))

    def set_values(self, idx, values):
        """Set the values of a variable.

        Parameters
        ----------
        idx : integer
            Index of the variable in the data array.
        values : numpy.ndarray
            2D array with the values to set.

        """
        self.data[:, :, idx] = values

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
        data = self.data[:, :, idx]
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
        self.data[i, j, idx] = value

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
