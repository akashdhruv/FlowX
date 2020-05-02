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

def apply_ibm(gridx, gridy, ivar, ibm_r):
    """Identify the points of the cylindrical solid body.

    Parameters
    ----------
    gridx : flowx.GridFaceX object
        The grid for the u velocity.
    gridy : flowx.GridFaceY object
        The grid for the v velocity.    
    ivar : string
        Name of the variable in the Grid structure.
    ibm_r : float
        radius of the cylinder
    
    """

    #subdomain grid facex
    x_ind_facex = numpy.zeros(2)
    listx = numpy.round(gridx.x.tolist(),4).tolist()
    x_ind_facex[0] = listx.index(-ibm_r)
    x_ind_facex[1] = listx.index(ibm_r)
    x_facex = gridx.x[int(x_ind_facex[0]):int(x_ind_facex[1]) + 1]

    y_ind_facex = numpy.zeros(2)
    listy = numpy.round(gridx.y.tolist(),4).tolist()
    y_ind_facex[0] = listy.index(-ibm_r - gridx.dy / 2)
    y_ind_facex[1] = listy.index(ibm_r + gridx.dy / 2)
    y_facex = gridx.y[int(y_ind_facex[0]):int(y_ind_facex[1]) + 1]

    #subdomain grid facey
    x_ind_facey = numpy.zeros(2)
    listx = numpy.round(gridy.x.tolist(),4).tolist()
    x_ind_facey[0] = listx.index(-ibm_r - gridx.dy / 2)
    x_ind_facey[1] = listx.index(ibm_r + gridx.dy / 2)
    x_facey = gridy.x[int(x_ind_facey[0]):int(x_ind_facey[1]) + 1]

    y_ind_facey = numpy.zeros(2)
    listy = numpy.round(gridy.y.tolist(),4).tolist()
    y_ind_facey[0] = listy.index(-ibm_r)
    y_ind_facey[1] = listy.index(ibm_r)
    y_facey = gridy.y[int(y_ind_facey[0]):int(y_ind_facey[1]) + 1]

    X_facex, Y_facex = numpy.meshgrid(x_facex, y_facex)
    X_facey, Y_facey = numpy.meshgrid(x_facey, y_facey)

    #distance from the centre of the cyrcle (xo = 0, yo = 0)
    d_facex = numpy.sqrt(X_facex**2 + Y_facex**2)
    d_facex = numpy.round(d_facex,2)

    d_facey = numpy.sqrt(X_facey**2 + Y_facey**2)
    d_facey = numpy.round(d_facey,2)

    #identification of points where d <= R
    temp = (d_facex <= ibm_r)
    indy_facex, indx_facex = numpy.nonzero(temp)

    temp = (d_facey <= ibm_r)
    indy_facey, indx_facey = numpy.nonzero(temp)

    #apply zero velocity on the identified points
    u  = gridx.get_values(ivar)
    v  = gridy.get_values(ivar)

    frac_u = numpy.ones(numpy.shape(d_facex))
    frac_v = numpy.ones(numpy.shape(d_facey))

    frac_u[indy_facex, indx_facex] = 0
    frac_v[indy_facey, indx_facey] = 0

    u[int(y_ind_facex[0]):int(y_ind_facex[1])+1, int(x_ind_facex[0]): int(x_ind_facex[1])+1] = (u[int(y_ind_facex[0]):int(y_ind_facex[1])+1, int(x_ind_facex[0]): 
                                int(x_ind_facex[1])+1] * frac_u)

    v[int(y_ind_facey[0]):int(y_ind_facey[1])+1, int(x_ind_facey[0]): int(x_ind_facey[1])+1] = (v[int(y_ind_facey[0]):int(y_ind_facey[1])+1, int(x_ind_facey[0]):
                                int(x_ind_facey[1])+1] * frac_v)
    return 

