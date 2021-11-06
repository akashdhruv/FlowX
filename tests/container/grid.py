"""Tests for `flowx/container/grid.py`."""

import numpy
import random
import unittest

import flowx.develop as flowx

class TestGrid(unittest.TestCase):
    """Unit-tests for the Grid class."""

    def setUp(self):
        """Set up the grid with 3 cell-centered variables."""
        self.num = 3
        self.varlist = ['ivar', 'asol', 'eror']
        self.nx, self.ny = 10, 20
        self.xblocks, self.yblocks = 1,1
        self.xmin, self.xmax = 0.0, 1.0
        self.ymin, self.ymax = -0.5, 0.5
        self.domain = flowx.domain.Domain(self.nx, self.ny, self.xblocks, self.yblocks,
                                          self.xmin, self.xmax, self.ymin, self.ymax, self.varlist)
        self.grid = self.domain[0]

    def test_getters_setters(self):
        """Test set/get values of a cell-centered variables."""
        varkey = self.varlist[0]
        values = numpy.ones((self.xblocks*self.yblocks,1,self.ny + 2,self.nx + 2))
        for block,index in zip(self.grid.blocklist,[*range(self.xblocks*self.yblocks)]):
            block[varkey] = values[index]
        values2 = self.grid[varkey]
        self.assertTrue(numpy.allclose(values, values2, atol=1e-12))

    def test_get_set(self):
        """Test set/get value of a cell-centered variables."""
        varkey = self.varlist[0]
        value = 0.123456
        i, j, k = random.randint(0, self.nx + 1), random.randint(0, self.ny + 1), 0
        blocktag = 0
        self.grid.blocklist[blocktag][varkey][k,j,i] = value
        value2 = self.grid[varkey][blocktag][k,j,i]
        self.assertEqual(value, value2)

    def test_get_cell_centered_coordinates(self):
        """Test gridline coordinates of a cell-centered grid."""
        dx = (self.xmax - self.xmin) / self.nx
        dy = (self.ymax - self.ymin) / self.ny
        x_test = numpy.linspace(self.xmin - dx / 2, self.xmax + dx / 2,
                                num=self.nx + 2)
        self.assertTrue(numpy.allclose(self.grid.blocklist[0].x, x_test, atol=1e-8))
        y_test = numpy.linspace(self.ymin - dy / 2, self.ymax + dy / 2,
                                num=self.ny + 2)
        self.assertTrue(numpy.allclose(self.grid.blocklist[0].y, y_test, atol=1e-8))

    def test_compute_error(self):
        """Test the error between two identical data arrays is zero."""
        ivar, asol, eror = self.varlist
        values = numpy.random.rand(1, self.ny + 2, self.nx + 2)
        self.grid.blocklist[0][ivar] = values
        self.grid.blocklist[0][asol] = values
        self.grid.compute_error(ivar, asol, eror)
        errors = self.grid.blocklist[0][eror]
        expected_errors = numpy.zeros((1, self.ny + 2, self.nx + 2))
        self.assertTrue(numpy.allclose(errors, expected_errors, atol=1e-12))

    def test_get_l2_norm(self):
        """Compare the L2-norm implementation with the NumPy implementation."""
        ivar = self.varlist[0]
        values = numpy.random.rand(1, self.ny + 2, self.nx + 2)
        self.grid.blocklist[0][ivar] = values
        l2_norm = self.grid.get_l2_norm(ivar)
        num = (self.nx + 2) * (self.ny + 2)
        expected_l2_norm = numpy.linalg.norm(values) / num
        self.assertTrue(abs(l2_norm - expected_l2_norm) < 1e-12)

if __name__ == '__main__':
    unittest.main()
