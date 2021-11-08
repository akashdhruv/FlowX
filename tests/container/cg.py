"""Tests for `flowx/poisson/cg.py`."""

import numpy
import random
import unittest

import flowx

class TestPoissonCG(unittest.TestCase):
    """Unit-tests for the Poisson CG solver."""

    def setUp(self):
        """Set up the grid and the variables of the Poisson system."""
        center_vars = ['ivar', 'rvar', 'asol', 'eror']
        nx, ny = 40, 40
        xmin, xmax = 0.0, 1.0
        ymin, ymax = -0.5, 0.5
        bc_type = {'ivar': 4 * ['dirichlet']}
        bc_val = {'ivar': 4 * [0.0]}
        self.grid = flowx.domain.Grid("cell-centered", center_vars, nx, ny, xmin, xmax, ymin, ymax,
                                       user_bc_type=bc_type, user_bc_val=bc_val)

        self._set_analytical('asol')
        self._set_rhs('rvar')

    def _set_analytical(self, var_name):
        """Private method to set the analytical solution."""
        X, Y = numpy.meshgrid(self.grid.x, self.grid.y)
        Lx = self.grid.xmax - self.grid.xmin
        Ly = self.grid.ymax - self.grid.ymin
        self.grid[var_name] = numpy.sin(numpy.pi * X / Lx) * numpy.cos(numpy.pi * Y / Ly)

    def _set_rhs(self, var_name):
        """Private method to set the right-hand side of the system."""
        X, Y = numpy.meshgrid(self.grid.x, self.grid.y)
        Lx = self.grid.xmax - self.grid.xmin
        Ly = self.grid.ymax - self.grid.ymin
        self.grid[var_name] = (-((numpy.pi / Lx)**2 + (numpy.pi / Ly)**2) *
                                  numpy.sin(numpy.pi * X / Lx) *
                             numpy.cos(numpy.pi * Y / Ly))

    def test_number_of_iterations(self):
        """Test the solver reaches the maximum number of iterations."""

        maxiter, tol = 0, 1e-12
        poisson_info = dict(poisson_solver='cg', maxiter=maxiter, tol=tol)
        poisson_vars = ['ivar', 'rvar']
        self.poisson = flowx.poisson.Poisson(self.grid, poisson_vars, poisson_info)
        ites, _ = self.poisson.solve()
        self.assertEqual(ites, maxiter)

    def test_residual(self):
        """Test the solver convergence."""
        maxiter, tol = 3000, 1e-6
        poisson_info = dict(poisson_solver='cg', maxiter=maxiter, tol=tol)
        poisson_vars = ['ivar', 'rvar']
        self.poisson = flowx.poisson.Poisson(self.grid, poisson_vars, poisson_info)
        ites, res = self.poisson.solve()
        self.assertTrue(res <= tol)
        self.assertTrue(ites < maxiter)

if __name__ == '__main__':
    unittest.main()
