"""Tests for `flowx/container/grid.py`."""

import numpy
import random
import unittest

import flowx

class TestGrid(unittest.TestCase):
    """Unit-tests for the Grid class."""

    def setUp(self):
        """Set up the grid with 3 cell-centered variables."""
        self.num = 3
        self.var_names = ['ivar', 'asol', 'eror']
        self.nx, self.ny = 10, 20
        self.xmin, self.xmax = 0.0, 1.0
        self.ymin, self.ymax = -0.5, 0.5
        self.domain = flowx.domain.Domain(self.nx, self.ny,
                                          self.xmin, self.xmax,
                                          self.ymin, self.ymax, self.var_names)

        self.grid = self.domain[0]

    def test_init(self):
        """Test initialization of attributes."""
        names = ['num', 'nx', 'ny', 'xmin', 'xmax', 'ymin', 'ymax']
        for name in names:
            self.assertEqual(getattr(self.grid, name), getattr(self, name))
        size = self.num * (self.nx + 2) * (self. ny + 2)
        self.assertEqual(self.grid.num * (self.grid.nx + 2) * (self.grid.ny + 2), size)

    def test_set_get_values(self):
        """Test set/get values of a cell-centered variables."""
        var_name = self.var_names[0]
        values = numpy.ones((self.nx + 2, self.ny + 2))
        self.grid.set_values(var_name, values)
        values2 = self.grid.get_values(var_name)
        self.assertTrue(numpy.allclose(values, values2, atol=1e-12))

    def test_set_get_value(self):
        """Test set/get value of a cell-centered variables."""
        var_name = self.var_names[0]
        value = 0.123456
        i, j = random.randint(0, self.nx + 1), random.randint(0, self.ny + 1)
        self.grid.set_value(var_name, i, j, value)
        value2 = self.grid.get_value(var_name, i, j)
        self.assertEqual(value, value2)

    def test_get_cell_centered_coordinates(self):
        """Test gridline coordinates of a cell-centered grid."""
        dx = (self.xmax - self.xmin) / self.nx
        dy = (self.ymax - self.ymin) / self.ny
        x_test = numpy.linspace(self.xmin - dx / 2, self.xmax + dx / 2,
                                num=self.nx + 2)
        self.assertTrue(numpy.allclose(self.grid.x, x_test, atol=1e-8))
        y_test = numpy.linspace(self.ymin - dy / 2, self.ymax + dy / 2,
                                num=self.ny + 2)
        self.assertTrue(numpy.allclose(self.grid.y, y_test, atol=1e-8))

    def test_get_error(self):
        """Test the error between two identical data arrays is zero."""
        ivar, asol, eror = self.var_names
        values = numpy.random.rand(self.nx + 2, self.ny + 2)
        self.grid.set_values(ivar, values)
        self.grid.set_values(asol, values)
        self.grid.get_error(ivar, asol, eror)
        errors = self.grid.get_values(eror)
        expected_errors = numpy.zeros((self.nx + 2, self.ny + 2))
        self.assertTrue(numpy.allclose(errors, expected_errors, atol=1e-12))

    def test_get_l2_norm(self):
        """Compare the L2-norm implementation with the NumPy implementation."""
        ivar = self.var_names[0]
        values = numpy.random.rand(self.nx + 2, self.ny + 2)
        self.grid.set_values(ivar, values)
        l2_norm = self.grid.get_l2_norm(ivar)
        num = (self.nx + 2) * (self.ny + 2)
        expected_l2_norm = numpy.linalg.norm(values) / num
        self.assertTrue(abs(l2_norm - expected_l2_norm) < 1e-12)


if __name__ == '__main__':
    unittest.main()
