import numpy
import time

from . import _interface


def advect(gridc, gridx, gridy, scalars, ibmf, ibmx, ibmy, velc, options):

    """ """

    nx, ny = gridc.nx, gridc.ny
    dx, dy = gridc.dx, gridc.dy
    dt = scalars.dt
    lset_iter = options["lset_redistance"]
    extrap_iter = options["extrap_solid"]

    u = gridx[velc][0, 0, :, :].transpose()
    v = gridy[velc][0, 0, :, :].transpose()

    phi = gridc[ibmf][0, 0, :, :].transpose()
    lmx = gridc[ibmx][0, 0, :, :].transpose()
    lmy = gridc[ibmy][0, 0, :, :].transpose()

    adfx = numpy.zeros_like(phi)
    adfy = numpy.zeros_like(phi)
    ddsn = numpy.zeros_like(phi)

    # --------------Advect dynamic X-Y grid---------------
    _interface.advect_dynamic_grid(phi, lmx, u, v, dt, dx, dy, nx + 2, ny + 2)
    _interface.advect_dynamic_grid(phi, lmy, u, v, dt, dx, dy, nx + 2, ny + 2)

    # ---------Find normal vectors---------------------
    _interface.normal_vector_solid(phi, adfx, adfy, dx, dy, nx + 2, ny + 2)

    # ---------Extrapolate X grid----------------------
    _interface.directional_derivative(
        phi, lmx, ddsn, adfx, adfy, dx, dy, nx + 2, ny + 2
    )

    for _iter in range(extrap_iter):
        _interface.constant_extrapolation(phi, ddsn, adfx, adfy, dx, dy, nx + 2, ny + 2)

    for _iter in range(extrap_iter):
        _interface.linear_extrapolation(
            phi, lmx, ddsn, adfx, adfy, dx, dy, nx + 2, ny + 2
        )
        gridc.fill_guard_cells(ibmx)

    # ---------Extrapolate Y grid----------------------
    _interface.directional_derivative(
        phi, lmy, ddsn, adfx, adfy, dx, dy, nx + 2, ny + 2
    )

    for _iter in range(extrap_iter):
        _interface.constant_extrapolation(phi, ddsn, adfx, adfy, dx, dy, nx + 2, ny + 2)

    for _iter in range(extrap_iter):
        _interface.linear_extrapolation(
            phi, lmy, ddsn, adfx, adfy, dx, dy, nx + 2, ny + 2
        )
        gridc.fill_guard_cells(ibmy)

    # --------Advect solid interface-------------------
    _interface.advect_solid(phi, u, v, dt, dx, dy, nx + 2, ny + 2)
    gridc.fill_guard_cells(ibmf)

    # --------Redistance solid interface---------------
    phi_old = numpy.copy(phi)
    lsDT = numpy.sqrt(dx**2 + dy**2) / 2.0

    for _iter in range(lset_iter):
        _interface.redistance_solid(phi, phi_old, lsDT, dx, dy, nx + 2, ny + 2)
        gridc.fill_guard_cells(ibmf)

    return
