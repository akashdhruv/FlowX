import numpy
from scipy import special


def solid_props(s, xmus, mu_solid, dx, dy, nx, ny):

    eps = 1.0e-10

    phi = -numpy.copy(s)

    psi = (1.0 + special.erf(-phi / (2.0 * dx))) / 2.0

    xmus[:, :] = psi * (mu_solid - 0.0) + 0.0

    return
