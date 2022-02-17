import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg.dsolve import linsolve
from scipy.sparse import spdiags, csr_matrix
from scipy.sparse import linalg as sla


def build_sparse_matrix(grid, ivar):

    """Solve the Poisson system using a direct solver from the scipy library.

     Arguments
     ---------
     grid : Grid object
         Grid containing data.

    Returns
     -------
     matrix : CSR format matrix

    """
    nx, ny = grid.nx, grid.ny
    dx, dy = grid.dx, grid.dy

    matrix_length = nx * ny

    matrix = sps.lil_matrix((matrix_length, matrix_length), dtype=np.float64)

    counter = 0

    bc_type = grid.bc_type[ivar]
    coeff_add = [None] * 4

    for i in range(len(bc_type)):
        if bc_type[i] == "neumann":
            coeff_add[i] = 1.0
        if bc_type[i] == "dirichlet":
            coeff_add[i] = -1.0

    coeff_start = -4.0

    for j in range(1, ny + 1):
        for i in range(1, nx + 1):

            coeff = coeff_start / dx**2

            if j > 1:
                matrix[counter, counter - nx] = 1.0 / (dx**2)
            else:
                coeff = coeff + coeff_add[2] / (dx**2)

            if j < ny:
                matrix[counter, counter + nx] = 1.0 / (dx**2)
            else:
                coeff = coeff + coeff_add[3] / (dx**2)

            if i > 1:
                matrix[counter, counter - 1] = 1.0 / (dx**2)
            else:
                coeff = coeff + coeff_add[0] / (dx**2)

            if i < nx:
                matrix[counter, counter + 1] = 1.0 / (dx**2)
            else:
                coeff = coeff + coeff_add[1] / (dx**2)

            matrix[counter, counter] = coeff
            counter = counter + 1

    matrix = matrix.tocsr()
    lu = sla.splu(matrix.tocsc())

    return lu, matrix
