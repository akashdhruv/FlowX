import numpy
from numba import jit


def solid_ustar(ustr, vstr, xms, Tau1, Tau2, Tau3, Tau4, Re_s, dt, dx, dy, nx, ny):

    _jit_solid_ustar(ustr, vstr, xms, Tau1, Tau2, Tau3, Tau4, Re_s, dt, dx, dy, nx, ny)

    return


@jit(nopython=True)
def _jit_solid_ustar(ustr, vstr, xms, Tau1, Tau2, Tau3, Tau4, Re_s, dt, dx, dy, nx, ny):

    for i in range(1, nx - 2):
        for j in range(1, ny - 1):

            txplus = xms[i + 1, j] * Tau1[i + 1, j]

            txmins = xms[i, j] * Tau1[i, j]

            typlus = (
                (xms[i + 1, j + 1] + xms[i, j + 1])
                / 2.0
                * (Tau2[i + 1, j + 1] + Tau2[i, j + 1])
                / 2.0
            )

            tymins = (
                (xms[i + 1, j - 1] + xms[i, j - 1])
                / 2.0
                * (Tau2[i + 1, j - 1] + Tau2[i, j - 1])
                / 2.0
            )

            # ---Update ustar by adding elastic solid force term to solid region

            ustrB = (txplus - txmins) / (dx * Re_s) + (typlus - tymins) / (
                2 * dy * Re_s
            )

            ustr[i, j] = ustr[i, j] + ustrB * dt

    for i in range(1, nx - 1):
        for j in range(1, ny - 2):

            txplus = (
                (xms[i + 1, j - 1] + xms[i + 1, j])
                / 2.0
                * (Tau3[i + 1, j - 1] + Tau3[i + 1, j])
                / 2.0
            )

            txmins = (
                (xms[i - 1, j - 1] + xms[i - 1, j])
                / 2.0
                * (Tau3[i - 1, j - 1] + Tau3[i - 1, j])
                / 2.0
            )

            typlus = xms[i, j + 1] * Tau4[i, j + 1]

            tymins = xms[i, j] * Tau4[i, j]

            # ---Update vstar by adding elastic solid force term to solid region

            vstrB = (txplus - txmins) / (2 * dx * Re_s) + (typlus - tymins) / (
                dy * Re_s
            )

            vstr[i, j] = vstr[i, j] + vstrB * dt

    return
