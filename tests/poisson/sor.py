"""Tests for `mae6225/poisson/sor.py`."""

import numpy
import random
import unittest

import mae6225


class TestPoissonSOR(unittest.TestCase):
    """Unit-tests for the Poisson SOR solver."""

    def setUp(self):
        """Set up the grid and the variables of the Poisson system."""
        center_vars = ['ivar', 'rvar', 'asol', 'eror']
        nx, ny = 40, 40
        xmin, xmax = 0.0, 1.0
        ymin, ymax = -0.5, 0.5
        bc_type = {'ivar': 4 * ['dirichlet']}
        bc_val = {'ivar': 4 * [0.0]}
        self.grid = mae6225.Grid('cell-centered',
                                 center_vars,
                                 nx, ny,
                                 xmin, xmax, ymin, ymax,
                                 user_bc_type=bc_type, user_bc_val=bc_val)
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
        ites, _ = mae6225.poisson.solve_sor(self.grid, 'ivar', 'rvar',
                                            maxiter=maxiter, tol=tol)
        self.assertEqual(ites, maxiter)
        maxiter, tol = 100, 1e-12
        ites, _ = mae6225.poisson.solve_sor(self.grid, 'ivar', 'rvar',
                                            maxiter=maxiter, tol=tol)
        self.assertEqual(ites, maxiter)

    def test_residual(self):
        """Test the solver convergence."""
        maxiter, tol = 3000, 1e-6
        ites, res = mae6225.poisson.solve_sor(self.grid, 'ivar', 'rvar',
                                              maxiter=maxiter, tol=tol)
        self.assertTrue(res <= tol)
        self.assertTrue(ites < maxiter)

    def test_error(self):
        """Test the solver results as we decrease the exit tolerance."""
        maxiter, tol = 3000, 1e-3
        ites1, res1 = mae6225.poisson.solve_sor(self.grid, 'ivar', 'rvar',
                                                maxiter=maxiter, tol=tol)
        self.grid.get_error('eror', 'ivar', 'asol')
        error1 = numpy.max(numpy.abs(self.grid.get_values('eror')))
        # Reset numerical solutio for second solve
        self.grid.set_values('ivar', numpy.zeros((self.grid.nx + 2,
                                                  self.grid.ny + 2)))
        maxiter, tol = 3000, 1e-6
        ites2, res2 = mae6225.poisson.solve_sor(self.grid, 'ivar', 'rvar',
                                                maxiter=maxiter, tol=tol)
        self.grid.get_error('eror', 'ivar', 'asol')
        error2 = numpy.max(numpy.abs(self.grid.get_values('eror')))
        self.assertTrue(ites1 <= ites2)
        self.assertTrue(res1 >= res2)
        self.assertTrue(error1 >= error2)


if __name__ == '__main__':
    unittest.main()
