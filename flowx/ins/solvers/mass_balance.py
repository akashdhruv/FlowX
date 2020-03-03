import numpy

def get_qin(grid, ivar):
    """Compute and return the mass getting in the domain.

    Arguments
    ---------
    grid : flowx.Grid object
        Grid containing data.
    ivar : string
        Name of the velocity variable on the grid.
    bctype : dictionary
        Type of boundary conditions for the variable `ivar`.

    Returns
    -------
    Qin : float
        Mass getting in the domain.

    """
    vel = grid.get_values(ivar)
    dx, dy = grid.dx, grid.dy

    bc_type = grid.bc_type[ivar]

    Qin = 0.0

    if grid.type_ == 'x-face':
        if bc_type[0] is not 'outflow' and bc_type[0] is not 'neumann':
            Qin += numpy.sum(vel[0, 1:-1]) * dy
        if bc_type[1] is not 'outflow' and bc_type[1] is not 'neumann':
            Qin -= numpy.sum(vel[-1, 1:-1]) * dy
    elif grid.type_ == 'y-face':
        if bc_type[2] is not 'outflow' and bc_type[2] is not 'neumann':
            Qin += numpy.sum(vel[1:-1, 0]) * dx
        if bc_type[3] is not 'outflow' and bc_type[3] is not 'neumann':
            Qin -= numpy.sum(vel[1:-1, -1]) * dx

    return Qin

def get_qout(grid, ivar):
    """Compute and return the mass getting out the domain.

    Arguments
    ---------
    grid : flowx.Grid object
        Grid containing data.
    ivar : string
        Name of the velocity variable on the grid.
    bctype : dictionary
        Type of boundary conditions for the variable `ivar`.

    Returns
    -------
    Qout : float
        Mass getting out the domain.

    """
    vel = grid.get_values(ivar)
    dx, dy = grid.dx, grid.dy

    bc_type = grid.bc_type[ivar]

    Qout = 0.0

    if grid.type_ == 'x-face':
        if bc_type[0] is 'outflow' or bc_type[0] is 'neumann':
            Qout -= numpy.sum(vel[0, 1:-1]) * dy
        if bc_type[1] is 'outflow' or bc_type[1] is 'neumann':
            Qout += numpy.sum(vel[-1, 1:-1]) * dy
    elif grid.type_ == 'y-face':
        if bc_type[2] is 'outflow' or bc_type[2] is 'neumann':
            Qout -= numpy.sum(vel[1:-1, 0]) * dx
        if bc_type[3] is 'outflow' or bc_type[3] is 'neumann':
            Qout += numpy.sum(vel[1:-1, -1]) * dx

    return Qout

def rescale_velocity(grid, ivar, Qin, Qout):
    """Rescale velocity.

    Arguments
    ---------
    grid : flowx.Grid object
        Grid containing data.
    ivar : string
        Name of the velocity variable on the grid.
    bctype : dictionary
        Type of boundary conditions for the variable `ivar`.
    Qin : float
        Mass in.
    Qout : float
        Mass out.

    """
    vel = grid.get_values(ivar)

    bc_type = grid.bc_type[ivar]

    Qinout = 1.0
    if Qout > 0.0:
        Qinout = Qin/Qout

    if grid.type_ == 'x-face':
        if bc_type[0] is 'outflow' or bc_type[0] is 'neumann':
            vel[0, 1:-1] *= Qinout
        if bc_type[1] is 'outflow' or bc_type[1] is 'neumann':
            vel[-1, 1:-1] *= Qinout

    if grid.type_ == 'y-face':
        if bc_type[2] is 'outflow' or bc_type[2] is 'neumann':
            vel[1:-1, 0] *= Qinout
        if bc_type[3] is 'outflow' or bc_type[3] is 'neumann':
            vel[1:-1, -1] *= Qinout

    return

def get_convvel(grid, ivar):
    """Get convective outflow velocity.

    Arguments
    ---------
    grid : flowx.Grid object
        Grid containing data.
    ivar : string
        Name of the velocity variable on the grid.

    Returns
    -------
    convvel : float
        Variable containing outflow velocity.

    """
    vel = grid.get_values(ivar)

    bc_type = grid.bc_type[ivar]

    convvel = [0.0, 0.0, 0.0, 0.0]

    if grid.type_ == 'x-face':
        if bc_type[0] is 'outflow':
            convvel[0] = numpy.mean(vel[0, :])
        if bc_type[1] is 'outflow':
            convvel[1] = numpy.mean(vel[-1, :])

    if grid.type_ == 'y-face':
        if bc_type[2] is 'outflow':
            convvel[2] = numpy.mean(vel[:, 0])
        if bc_type[3] is 'outflow':
            convvel[3] = numpy.mean(vel[:, -1])

    return convvel

def update_outflow_bc(grid, ivar, dt, convvel=None):
    """Update the value of the velocity at the right boundary.

    The function uses a linear convective equation in the x-direction
    where the convective velocity is defined as the mean x-velocity
    along the right boundary.

    Parameters
    ----------
    grid : flowx.GridFaceX object
        The grid for the velocity.
    ivar : string
        Name of the variable in the Grid structure.
    dt : float
       Time-step size.
    convvel : float (optional)
        Convective velocity;
        default: None (will compute the convective velocity).

    """
    vel = grid.get_values(ivar)
    dx,dy = grid.dx, grid.dy

    bc_type = grid.bc_type[ivar]
    bc_val = grid.bc_val[ivar]

    grid_type = grid.type_

    if convvel is None:
        convvel = get_convvel(grid, ivar)

    if grid.type_ == 'x-face':
        if bc_type[0] is 'outflow' or bc_type[0] is 'neumann':
            bc_val[0] = vel[0, :] - convvel[0] * dt * (vel[1, :] - vel[0, :]) / dx

        if bc_type[1] is 'outflow' or bc_type[1] is 'neumann':
            bc_val[1] = vel[-1, :] - convvel[1] * dt * (vel[-1, :] - vel[-2, :]) / dx

    if grid.type_ == 'y-face':
        if bc_type[2] is 'outflow' or bc_type[2] is 'neumann':
            bc_val[2] = vel[:, 0] - convvel[2] * dt * (vel[:, 1] - vel[:, 0]) / dy

        if bc_type[3] is 'outflow' or bc_type[3] is 'neumann':
            bc_val[3] = vel[:, -1] - convvel[3] * dt * (vel[:, -1] - vel[:, -2]) / dy

    grid.update_bc_val({ivar: bc_val})

    return
