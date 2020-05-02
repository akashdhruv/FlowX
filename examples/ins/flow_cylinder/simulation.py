"""User-defined module for simulation."""

import numpy


def set_initial_velocity(gridc, gridx, gridy, ivar, pres):
    """Set the initial velocity field.

    The x- and y-components of the velocity are set to 1.0 and 0.0,
    respectively.

    Arguments
    ---------
    gridx : flowx.Grid object
        Grid containing x-face data.
    gridy : flowx.Grid object
        Grid containing y-face data.
    ivar : string
        Name of the velocity variable on the grid.

    """

    u = gridx.get_values(ivar)
    v = gridy.get_values(ivar)
    p = gridc.get_values(pres)

    u[:, :] = 1.0
    v[:, :] = 0.0
    p[:, :] = 0.0

    return


def get_qin(grid, ivar, bc_type):
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

    Qin = 0.0

    if grid.type_ is 'x-face':
        if bc_type[0] is not 'outflow':
            Qin += numpy.sum(vel[0, 1:-1]) * dy
        if bc_type[1] is not 'outflow':
            Qin -= numpy.sum(vel[-1, 1:-1]) * dy
    elif grid.type_ is 'y-face':
        if bc_type[2] is not 'outflow':
            Qin += numpy.sum(vel[1:-1, 0]) * dx
        if bc_type[3] is not 'outflow':
            Qin -= numpy.sum(vel[1:-1, -1]) * dx

    return Qin


def get_qout(grid, ivar, bc_type):
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

    Qout = 0.0

    if grid.type_ is 'x-face':
        if bc_type[0] is 'outflow':
            Qout -= numpy.sum(vel[0, 1:-1]) * dy
        if bc_type[1] is 'outflow':
            Qout += numpy.sum(vel[-1, 1:-1]) * dy
    elif grid.type_ is 'y-face':
        if bc_type[2] is 'outflow':
            Qout -= numpy.sum(vel[1:-1, 0]) * dx
        if bc_type[3] is 'outflow':
            Qout += numpy.sum(vel[1:-1, -1]) * dx

    return Qout


def rescale_velocity(grid, ivar, bc_type, Qin, Qout):
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

    Qinout = 1.0
    if Qout > 0.0:
        Qinout = Qin/Qout

    if grid.type_ is 'x-face':
        if bc_type[0] is 'outflow':
            vel[0, 1:-1] *= Qinout
        if bc_type[1] is 'outflow':
            vel[-1, 1:-1] *= Qinout

    if grid.type_ is 'y-face':
        if bc_type[2] is 'outflow':
            vel[1:-1, 0] *= Qinout
        if bc_type[3] is 'outflow':
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

    convvel = numpy.mean(vel[-1, :])

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
    dx = grid.dx

    if convvel is None:
        convvel = get_convvel(grid, ivar)

    bc_val = grid.bc_val[ivar]
    bc_val[1] = vel[-1, :] - convvel * dt * (vel[-1, :] - vel[-2, :]) / dx
    grid.update_bc_val({ivar: bc_val})

    return

def ibm_vof(gridx, gridy, ibm_r, ivar):
    """Impose the effect of the solid body (cylinder) on the flow.

    Parameters
    ----------
    gridx : flowx.GridFaceX object
        The grid for the u velocity component.
    gridy : flowx.GridFaceY object
        The grid for the v velocity component.
    ibm_r : float
       The radius of the cylinder.
    ivar : string
        Name of the variable in the Grid structure.
    """

    #For the face_x grid (u velocity component)
    x_ind_face_x = numpy.zeros(2)
    list_x = numpy.round(gridx.x.tolist(),4).tolist()
    x_ind_face_x[0] = list_x.index(-ibm_r)
    x_ind_face_x[1] = list_x.index(ibm_r)
    x_face_x = gridx.x[int(x_ind_face_x[0]):int(x_ind_face_x[1]) + 1]

    y_ind_face_x = numpy.copy(x_ind_face_x)
    list_y = numpy.round(gridx.y.tolist(),4).tolist()
    y_ind_face_x[0] = list_y.index(-ibm_r - gridx.dy / 2)
    y_ind_face_x[1] = list_y.index(ibm_r + gridx.dy / 2)
    y_face_x = gridx.y[int(y_ind_face_x[0]):int(y_ind_face_x[1]) + 1]

    X_grid_face_x, Y_grid_face_x = numpy.meshgrid(x_face_x, y_face_x)
    d_face_x = numpy.sqrt(X_grid_face_x**2 + Y_grid_face_x**2)

    check = (d_face_x <= ibm_r)
    ind_y_face_x, ind_x_face_x = numpy.nonzero(check)
   
    u  = gridx.get_values(ivar)

    vf_u = numpy.ones(numpy.shape(d_face_x))
    vf_u[ind_y_face_x, ind_x_face_x] = 0

    u[int(y_ind_face_x[0]):int(y_ind_face_x[1])+1, int(x_ind_face_x[0]): int(x_ind_face_x[1])+1] = (u[int(y_ind_face_x[0]):
                           int(y_ind_face_x[1])+1, int(x_ind_face_x[0]): int(x_ind_face_x[1])+1] * vf_u)


    #For the face_y grid (v velocity component)
    x_ind_face_y = numpy.copy(x_ind_face_x)
    list_x = numpy.round(gridy.x.tolist(),4).tolist()
    x_ind_face_y[0] = list_x.index(-ibm_r - gridx.dy / 2)
    x_ind_face_y[1] = list_x.index(ibm_r + gridx.dy / 2)
    x_face_y = gridy.x[int(x_ind_face_y[0]):int(x_ind_face_y[1]) + 1]

    y_ind_face_y = numpy.copy(y_ind_face_x)
    list_y = numpy.round(gridy.y.tolist(),4).tolist()
    y_ind_face_y[0] = list_y.index(-ibm_r)
    y_ind_face_y[1] = list_y.index(ibm_r)
    y_face_y = gridy.y[int(y_ind_face_y[0]):int(y_ind_face_y[1]) + 1]

    X_grid_face_y, Y_grid_face_y = numpy.meshgrid(x_face_y, y_face_y)
    d_face_y = numpy.sqrt(X_grid_face_y**2 + Y_grid_face_y**2)

    check = (d_face_y <= ibm_r)
    ind_y_face_y, ind_x_face_y = numpy.nonzero(check)

    v  = gridy.get_values(ivar)
    
    vf_v = numpy.ones(numpy.shape(d_face_y))
    vf_v[ind_y_face_y, ind_x_face_y] = 0

    v[int(y_ind_face_y[0]):int(y_ind_face_y[1])+1, int(x_ind_face_y[0]): int(x_ind_face_y[1])+1] = (v[int(y_ind_face_y[0]):
                           int(y_ind_face_y[1])+1, int(x_ind_face_y[0]): int(x_ind_face_y[1])+1] * vf_v)


    return


