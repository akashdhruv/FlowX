"""Tests for `mae6225/poisson/jacobi.py`."""

import numpy
import random
import unittest

import mae6225


class TestPoissonJacobi(unittest.TestCase):
    def setUp(self):
        center_vars = ['ivar', 'rvar', 'asol', 'eror']
        nx, ny = 40, 40
        xmin, xmax = 0.0, 1.0
        ymin, ymax = -0.5, 0.5
        bc_type = {'ivar': 4 * ['dirichlet']}
        bc_val = {'ivar': 4 * [0.0]}
        self.grid = mae6225.Grid(center_vars,
                                 nx, ny,
                                 xmin, xmax, ymin, ymax,
                                 user_bc_type=bc_type, user_bc_val=bc_val)
        self._set_analytical('asol')
        self._set_rhs('rvar')

    def _set_analytical(self, var_name):
        x_center, y_center = self.grid.get_cell_centered_coordinates()
        X, Y = numpy.meshgrid(x_center, y_center)
        Lx = self.grid.xmax - self.grid.xmin
        Ly = self.grid.ymax - self.grid.ymin
        values = numpy.sin(numpy.pi * X / Lx) * numpy.cos(numpy.pi * Y / Ly)
        self.grid.set_values(var_name, values)

    def _set_rhs(self, var_name):
        x_center, y_center = self.grid.get_cell_centered_coordinates()
        X, Y = numpy.meshgrid(x_center, y_center)
        Lx = self.grid.xmax - self.grid.xmin
        Ly = self.grid.ymax - self.grid.ymin
        values = (-((numpy.pi / Lx)**2 + (numpy.pi / Ly)**2) *
                  numpy.sin(numpy.pi * X / Lx) *
                  numpy.cos(numpy.pi * Y / Ly))
        self.grid.set_values(var_name, values)

    def test_number_of_iterations(self):
        maxiter, tol = 0, 1e-12
        ites, _ = mae6225.poisson.solve_jacobi(self.grid, 'ivar', 'rvar',
                                               maxiter=maxiter, tol=tol)
        self.assertEqual(ites, maxiter)
        maxiter, tol = 100, 1e-12
        ites, _ = mae6225.poisson.solve_jacobi(self.grid, 'ivar', 'rvar',
                                               maxiter=maxiter, tol=tol)
        self.assertEqual(ites, maxiter)

    def test_residual(self):
        maxiter, tol = 3000, 1e-6
        ites, res = mae6225.poisson.solve_jacobi(self.grid, 'ivar', 'rvar',
                                                 maxiter=maxiter, tol=tol)
        self.assertTrue(res <= tol)
        self.assertTrue(ites < maxiter)

    def test_error(self):
        maxiter, tol = 3000, 1e-3
        ites1, res1 = mae6225.poisson.solve_jacobi(self.grid, 'ivar', 'rvar',
                                                   maxiter=maxiter, tol=tol)
        mae6225.poisson.get_error(self.grid, 'eror', 'ivar', 'asol')
        error1 = numpy.max(numpy.abs(self.grid.get_values('eror')))
        # Reset numerical solutio for second solve
        self.grid.set_values('ivar', numpy.zeros((self.grid.nx + 2,
                                                  self.grid.ny + 2)))
        maxiter, tol = 3000, 1e-6
        ites2, res2 = mae6225.poisson.solve_jacobi(self.grid, 'ivar', 'rvar',
                                                   maxiter=maxiter, tol=tol)
        mae6225.poisson.get_error(self.grid, 'eror', 'ivar', 'asol')
        error2 = numpy.max(numpy.abs(self.grid.get_values('eror')))
        self.assertTrue(ites1 <= ites2)
        self.assertTrue(res1 >= res2)
        self.assertTrue(error1 >= error2)


if __name__ == '__main__':
    unittest.main()
