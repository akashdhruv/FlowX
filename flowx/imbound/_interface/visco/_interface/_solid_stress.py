import numpy
from numba import jit


def solid_stress(sd, sX, sY, Tau1, Tau2, Tau3, Tau4, dx, dy, nx, ny):

    Tau1[:, :] = 0.0
    Tau2[:, :] = 0.0
    Tau3[:, :] = 0.0
    Tau4[:, :] = 0.0

    jit_solid_stress(sd, sX, sY, Tau1, Tau2, Tau3, Tau4, dx, dy, nx, ny)

    Tau1[0, :] = Tau1[1, :]
    Tau1[-1, :] = Tau1[-2, :]
    Tau1[:, 0] = Tau1[:, 1]
    Tau1[:, -1] = Tau1[:, -2]

    Tau2[0, :] = Tau2[1, :]
    Tau2[-1, :] = Tau2[-2, :]
    Tau2[:, 0] = Tau2[:, 1]
    Tau2[:, -1] = Tau2[:, -2]

    Tau3[0, :] = Tau3[1, :]
    Tau3[-1, :] = Tau3[-2, :]
    Tau3[:, 0] = Tau3[:, 1]
    Tau3[:, -1] = Tau3[:, -2]

    Tau4[0, :] = Tau4[1, :]
    Tau4[-1, :] = Tau4[-2, :]
    Tau4[:, 0] = Tau4[:, 1]
    Tau4[:, -1] = Tau4[:, -2]


@jit(nopython=True)
def jit_solid_stress(sd, sX, sY, Tau1, Tau2, Tau3, Tau4, dx, dy, nx, ny):

    for i in range(1, nx - 1):
        for j in range(1, ny - 1):

            A1 = 0.0
            A2 = 0.0
            A3 = 0.0
            A4 = 0.0

            AT1 = 0.0
            AT2 = 0.0
            AT3 = 0.0
            AT4 = 0.0

            A_inv1 = 1.0 / (2 * dx) * (sX[i + 1, j] - sX[i - 1, j])
            A_inv2 = 1.0 / (2 * dy) * (sX[i, j + 1] - sX[i, j - 1])
            A_inv3 = 1.0 / (2 * dx) * (sY[i + 1, j] - sY[i - 1, j])
            A_inv4 = 1.0 / (2 * dy) * (sY[i, j + 1] - sY[i, j - 1])

            if numpy.abs(A_inv1 * A_inv4 - A_inv2 * A_inv3) > 1e-12:

                A1 = 1.0 / (A_inv1 * A_inv4 - A_inv2 * A_inv3) * A_inv4
                A2 = 1.0 / (A_inv1 * A_inv4 - A_inv2 * A_inv3) * (-A_inv2)
                A3 = 1.0 / (A_inv1 * A_inv4 - A_inv2 * A_inv3) * (-A_inv3)
                A4 = 1.0 / (A_inv1 * A_inv4 - A_inv2 * A_inv3) * A_inv1

            AT1 = A1
            AT2 = A3
            AT3 = A2
            AT4 = A4

            Tau1[i, j] = A1 * AT1 + A2 * AT3 - 1
            Tau2[i, j] = A1 * AT2 + A2 * AT4
            Tau3[i, j] = A3 * AT1 + A4 * AT3
            Tau4[i, j] = A3 * AT2 + A4 * AT4 - 1

    return
