"""User-defined module for simulation."""

import numpy
import flowx

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

    u[:, :] =  1.0
    v[:, :] =  0.0
    p[:, :] =  0.0

    return

def main():

    # Define grid parameters
    nx, ny = 200, 400
    xmin, xmax = -2.5, 2.5
    ymin, ymax = -5.0, 5.0

    # Define cell-centered variable names
    center_vars   = ['pres', 'divv']
    face_vars     = ['velc', 'hvar', 'ibmf']
    ins_vars      = ['velc', 'hvar', 'divv', 'pres']
    poisson_vars  = ['pres', 'divv']
    imbound_vars  = ['ibmf', 'velc']

    scalar_info   = dict(tmax = 8, dt = 0.000625, Re = 100.0)

    simulation_info = dict(time_stepping = 'ab2', 
                           poisson_solver = 'serial_lu', 
                           maxiter = 2000,
                           tol = 1e-10,
                       with_ib = True)

    particle_info = [dict(file='sm_body.00001.h5', vel = [0.0, -1.0])]

    # Define boundary conditions for variable pressure and velocity [left, right, bottom, top]
    bc_type_center = dict(pres = ['neumann', 'neumann', 'neumann', 'neumann'])
    bc_val_center  = dict(pres = [0.0, 0.0, 0.0, 0.0])

    bc_type_facex = dict(velc = ['dirichlet', 'dirichlet', 'dirichlet', 'dirichlet'])
    bc_val_facex  = dict(velc = [0.0, 0.0, 0.0, 0.0])

    bc_type_facey = dict(velc = ['dirichlet', 'dirichlet', 'dirichlet', 'dirichlet'])
    bc_val_facey  = dict(velc =[0.0, 0.0, 0.0, 0.0])

    # Create the grid and data
    gridc, gridx, gridy, scalars, particles = flowx.serial.domain_main(nx, ny, xmin, xmax, ymin, ymax, \
                                                   center_vars, face_vars, scalar_info, particle_info, \
                                           bc_type_center=bc_type_center, bc_val_center=bc_val_center, \
                                               bc_type_facex=bc_type_facex, bc_val_facex=bc_val_facex, \
                                               bc_type_facey=bc_type_facey, bc_val_facey=bc_val_facey)

    domain_data_struct = [gridc, gridx, gridy, scalars, particles]

    poisson = flowx.poisson_main(gridc, poisson_vars, simulation_info)

    imbound = flowx.imbound_main(domain_data_struct, imbound_vars, simulation_info)

    ins = flowx.ins_main(poisson, imbound, domain_data_struct, ins_vars, simulation_info)

    while scalars.time <= scalars.tmax:
    
        imbound.map_to_grid()

        ins.advance()

        for particle in particles: particle.advance(scalars)

        scalars.advance()
    
        # Display stats
        if scalars.nstep % 1 == 0: flowx.io.display_stats(scalars)  
        
if __name__ == "__main__":
    main()
