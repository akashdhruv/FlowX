"""User defined module for simulation."""
  
import numpy

def set_initial_velocity(gridx, gridy, ivar):
    """Set initial condition for velocity.

    Arguments
    ---------   
    gridx : mae6225.Grid object
        Grid containing x-face data
    gridy : mae6225.Grid object
        Grid containing y-face data
    ivar : string
        Name of the velocity variable on the grid.
    """

    u = gridx.get_values(ivar)
    v = gridy.get_values(ivar)

    u[:,:]  =  1.0
    v[:,:]  =  0.0

    return

def get_qinqout(grid, ivar, bctype, flag):
    """Get Qin and Qout.

    Arguments
    ---------
    grid : mae6225.Grid object
        Grid containing data.
    ivar : string
        Name of the velocity variable on the grid.
    bctype : bc_type at different boundaries
    flag : logical
        Qin/Qout flag

    Returns
    --------
    Qaux : float
        Qinout value
    """

    vel = grid.get_values(ivar)

    dx,dy = grid.dx, grid.dy

    Qaux = 0
   
    if(flag == True):
  
        if(grid.type_ == 'x-face'):

            if(bctype[ivar][0] != 'outflow'):
                Qaux = Qaux + numpy.sum(vel[0,1:-1])*dy

            if(bctype[ivar][1] != 'outflow'):
                Qaux = Qaux - numpy.sum(vel[-1,1:-1])*dy
    
        elif(grid.type_ == 'y-face'):

            if(bctype[ivar][2] != 'outflow'):
                Qaux = Qaux + numpy.sum(vel[1:-1,0])*dx

            if(bctype[ivar][3] != 'outflow'):
                Qaux = Qaux - numpy.sum(vel[1:-1,-1])*dx

    elif(flag == False):
  
        if(grid.type_ == 'x-face'):

            if(bctype[ivar][0] == 'outflow'):
                Qaux = Qaux - numpy.sum(vel[0,1:-1])*dy

            if(bctype[ivar][1] == 'outflow'):
                Qaux = Qaux + numpy.sum(vel[-1,1:-1])*dy
    
        elif(grid.type_ == 'y-face'):

            if(bctype[ivar][2] == 'outflow'):
                Qaux = Qaux - numpy.sum(vel[1:-1,0])*dx

            if(bctype[ivar][3] == 'outflow'):
                Qaux = Qaux + numpy.sum(vel[1:-1,-1])*dx
 
    return Qaux

def rescale_velocity(grid, ivar, bctype, Qin, Qout):
    """Rescale velocity.

    Arguments
    ---------
    grid : mae6225.Grid object
        Grid containing data.
    ivar : string
        Name of the velocity variable on the grid.
    bctype : bc_type at different boundaries
    Qin : float
        Mass in
    Qout : float
        Mass out
    """

    vel = grid.get_values(ivar)

    Qinout = 1.0

    if(Qout > 1e-13):
        Qinout = Qin/Qout

    if(grid.type_ == 'x-face'):

        if(bctype[ivar][0] == 'outflow'):
            vel[0,1:-1] = vel[0,1:-1]*Qinout

        if(bctype[ivar][1] == 'outflow'):
            vel[-1,1:-1] = vel[-1,1:-1]*Qinout

    if(grid.type_ == 'y-face'):
    
        if(bctype[ivar][2] == 'outflow'):
            vel[1:-1,0] = vel[1:-1,0]*Qinout

        if(bctype[ivar][3] == 'outflow'):
            vel[1:-1,-1] = vel[1:-1,-1]*Qinout
  
    return

def get_convvel(grid, ivar):
    """get convective outflow velocity.

    Arguments
    ---------
    grid : mae6225.Grid object
        Grid containing data.
    ivar : string
        Name of the velocity variable on the grid.

    Returns
    --------
    convvel : float
          Variable containing outflow velocity
    """

    vel = grid.get_values(ivar)

    convvel = numpy.mean(vel[-1,:])

    return convvel

def update_outflow_bc(grid, ivar, convvel=0., dt=1., flg=False):
    """Update Dirichlet boundary values for the velocity components.

    Parameters
    ----------
    grid : mae6225.GridFaceX object
        The grid for the velocity.
    ivar : string
        Name of the velocity variable in the Grid structures.
    convvel : float
        Convective velocity   
    dt : float
       Time step
    flg : logical
      Flag for pred-corr
    """

    vel = grid.get_values(ivar)

    dx,dy = grid.dx, grid.dy

    if(flg):
        bc_val = vel[-1,:] - convvel*(dt/dx)*(vel[-1,:]-vel[-2,:]) 
    else:
        bc_val = vel[-1,:]

    grid.update_bc_val({ivar: [1.0,bc_val,0.0,0.0]})

    return

def ibm_tag_body(grid,ivar,ibm_x,ibm_y,ibm_r):
    """Tag immersed boundary

    Arguments
    ---------
    grid : mae6225.Grid object
        Grid containing data.
    ivar : string
        Name of the ibm tagging variable on the grid.
    """

    return

def ibm_apply_forcing(grid,ivar,ibmf,forc,ifac):
    """Apply immersed boundary forcing

    Arguments
    ---------
    grid : mae6225.Grid object
        Grid containing data.
    ivar : string
        Name of the velocity variable on the grid.
    ibmf : string
        Name of the tagging variable
    forc : string
        Name of the forcing variable
    """
    returnn
