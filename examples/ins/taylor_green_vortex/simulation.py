"""User defined module for simulation."""
  
import numpy
import flowx

def get_initial(gridc, gridx, gridy, ivar, pvar):
    """Compute and set the analytical solution.

    Arguments
    ---------
    gridc : flowx.Grid object
        Grid containing cell-center data.
    
    gridx : flowx.Grid object
        Grid containing x-face data

    girdy : flowx.Grid object
        Grid containing y-face data

    ivar : string
        Name of the velocity variable on the grid.

    pvar : string
        Name of the pressure variable on the grid.

    """

    Xcc, Ycc = numpy.meshgrid(gridc.x, gridc.y)
    Xfx, Yfx = numpy.meshgrid(gridx.x, gridx.y)
    Xfy, Yfy = numpy.meshgrid(gridy.x, gridy.y)

    Xcc = Xcc.transpose()
    Ycc = Ycc.transpose()

    Xfx = Xfx.transpose()
    Yfx = Yfx.transpose()

    Xfy = Xfy.transpose()
    Yfy = Yfy.transpose()
 
    u = gridx.get_values(ivar)
    v = gridy.get_values(ivar)
    p = gridc.get_values(pvar)

    u[:,:] =  -numpy.cos(Xfx)*numpy.sin(Yfx)
    v[:,:] =   numpy.sin(Xfy)*numpy.cos(Yfy)
    p[:,:] = -(numpy.cos(2*Xcc)+numpy.sin(2*Ycc))/4

    return

def get_analytical(gridc, gridx, gridy, asol, ifac):
    """Compute and set the analytical solution.

    Arguments
    ---------
    gridc : flowx.Grid object
        Grid containing cell-center data.
    
    gridx : flowx.Grid object
        Grid containing x-face data

    girdy : flowx.Grid object
        Grid containing y-face data

    asol : string
        Name of the variable on the grid.

    ifac : float
       Time step

    """

    Xcc, Ycc = numpy.meshgrid(gridc.x, gridc.y)
    Xfx, Yfx = numpy.meshgrid(gridx.x, gridx.y)
    Xfy, Yfy = numpy.meshgrid(gridy.x, gridy.y)

    Xcc = Xcc.transpose()
    Ycc = Ycc.transpose()

    Xfx = Xfx.transpose()
    Yfx = Yfx.transpose()

    Xfy = Xfy.transpose()
    Yfy = Yfy.transpose()

    u = gridx.get_values(asol)
    v = gridy.get_values(asol)
    p = gridc.get_values(asol)

    u[:,:] =  -numpy.exp(-2*ifac)*numpy.cos(Xfx)*numpy.sin(Yfx)
    v[:,:] =   numpy.exp(-2*ifac)*numpy.sin(Xfy)*numpy.cos(Yfy)
    p[:,:] =  -numpy.exp(-4*ifac)*(numpy.cos(2*Xcc)+numpy.sin(2*Ycc))/4

    return

def update_bc_val(gridx, gridy, ivar, t):
    """Update Dirichlet boundary values for the velocity components.

    Parameters
    ----------
    gridx : flowx.GridFaceX object
        The grid for the x-component of the velocity.
    gridy : flowx.GridFaceY object
        The grid for the y-component of the velocity.
    ivar : string
        Name of the velocity variable in the Grid structures.
    t : float
        Time.
    """
    coeff = numpy.exp(-2 * t)
    bc_val_u = [-coeff * numpy.cos(gridx.xmin) * numpy.sin(gridx.y),
                -coeff * numpy.cos(gridx.xmax) * numpy.sin(gridx.y),
                -coeff * numpy.cos(gridx.x) * numpy.sin(gridx.ymin),
                -coeff * numpy.cos(gridx.x) * numpy.sin(gridx.ymax)]
    bc_val_v = [coeff * numpy.sin(gridy.xmin) * numpy.cos(gridy.y),
                coeff * numpy.sin(gridy.xmax) * numpy.cos(gridy.y),
                coeff * numpy.sin(gridy.x) * numpy.cos(gridy.ymin),
                coeff * numpy.sin(gridy.x) * numpy.cos(gridy.ymax)]
    gridx.update_bc_val({ivar: bc_val_u})
    gridy.update_bc_val({ivar: bc_val_v})
    return

def test():
    # Define grid parameters
    nx, ny = 40, 40

    xmin, xmax = 0.0, 2.0*numpy.pi
    ymin, ymax = 0.0, 2.0*numpy.pi

    # Define cell-centered variable names
    center_vars   = ['pres', 'delp', 'divv', 'asol', 'eror']
    face_vars     = ['velc', 'hvar', 'asol', 'eror']

    ins_vars      = ['velc', 'hvar', 'divv', 'pres', 'delp']
    poisson_vars  = ['delp', 'divv']

    scalar_info     = dict(tmax =  2, dt = 0.001, Re = 1.0)
    simulation_info = dict(time_stepping = 'euler', 
                       poisson_solver = 'direct', 
                       maxiter = 3000, 
                       pressure_correct = True)

    # Define boundary conditions for variable pressure and velocity [left, right, bottom, top]
    bc_type_center = dict(delp = ['neumann', 'neumann', 'neumann', 'neumann'])
    bc_val_center  = dict(delp = [0.0, 0.0, 0.0, 0.0])

    bc_type_facex = dict(velc = ['dirichlet', 'dirichlet', 'dirichlet', 'dirichlet'])
    bc_val_facex  = dict(velc = [0.0, 0.0, 0.0, 0.0])

    bc_type_facey = dict(velc = ['dirichlet', 'dirichlet', 'dirichlet', 'dirichlet'])
    bc_val_facey  = dict(velc = [0.0, 0.0, 0.0, 0.0])

    gridc, gridx, gridy, scalars, particles = flowx.domain.Domain(nx, ny, xmin, xmax, ymin, ymax,
                                              center_vars, face_vars, scalar_info,
                                              bc_type_center=bc_type_center, bc_val_center=bc_val_center,
                                              bc_type_facex=bc_type_facex, bc_val_facex=bc_val_facex,
                                              bc_type_facey=bc_type_facey, bc_val_facey=bc_val_facey)

    domain_data_struct = [gridc, gridx, gridy, scalars, particles]

    poisson = flowx.poisson.Poisson(gridc, poisson_vars, simulation_info)

    imbound = flowx.imbound.ImBound()

    ins = flowx.ins.IncompNS(poisson, imbound, domain_data_struct, ins_vars, simulation_info)

    update_bc_val(gridx, gridy, 'velc', scalars.to)
    get_initial(gridc, gridx, gridy, 'velc', 'pres')

    while scalars.time <= scalars.tmax:
    
        # Update the time-dependent boundary condition value
        update_bc_val(gridx, gridy, 'velc', scalars.time)
    
        ins.advance()
    
        # Display stats
        if scalars.nstep % 10 == 0:
            flowx.io.display_stats(scalars)   
    
        scalars.advance()

    get_analytical(gridc, gridx, gridy, 'asol', scalars.time)


    gridx.get_error('eror','velc','asol')
    gridy.get_error('eror','velc','asol')

    #print(numpy.max(gridx.get_values('eror')),numpy.min(gridx.get_values('eror')))
    #print(numpy.max(gridy.get_values('eror')),numpy.min(gridy.get_values('eror')))

    maxdivv,mindivv = numpy.max(gridc.get_values('divv')),numpy.min(gridc.get_values('divv'))

    if (abs(maxdivv) <= 1e-11 and abs(mindivv) <= 1e-11 and maxdivv*mindivv < 0.):
        print('Divergence is within tolerance')
    else:
        raise ValueError('Divergence is not within tolerance')

if __name__ == "__main__":
    test()
