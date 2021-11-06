"""Tests for `flowx/poisson/cg.py`."""

import numpy
import random
import unittest

import flowx.archive as flowx

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
        self.domain = flowx.domain.Domain(nx, ny,
                                          xmin, xmax, ymin, ymax, center_vars,
                                          bc_type_center=bc_type, bc_val_center=bc_val)
        self.grid = self.domain[0]

        self._set_analytical('asol')
        self._set_rhs('rvar')

    def _set_analytical(self, var_name):
        """Private method to set the analytical solution."""
        X, Y = numpy.meshgrid(self.grid.x, self.grid.y)
        Lx = self.grid.xmax - self.grid.xmin
        Ly = self.grid.ymax - self.grid.ymin
        values = numpy.sin(numpy.pi * X / Lx) * numpy.cos(numpy.pi * Y / Ly)
        self.grid.set_values(var_name, values)

    def _set_rhs(self, var_name):
        """Private method to set the right-hand side of the system."""
        X, Y = numpy.meshgrid(self.grid.x, self.grid.y)
        Lx = self.grid.xmax - self.grid.xmin
        Ly = self.grid.ymax - self.grid.ymin
        values = (-((numpy.pi / Lx)**2 + (numpy.pi / Ly)**2) *
                  numpy.sin(numpy.pi * X / Lx) *
                  numpy.cos(numpy.pi * Y / Ly))
        self.grid.set_values(var_name, values)

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
