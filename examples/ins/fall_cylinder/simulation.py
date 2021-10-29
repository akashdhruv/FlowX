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
    nx, ny = 40, 80
    xmin, xmax = -2.5, 2.5
    ymin, ymax = -5.0, 5.0

    # Define cell-centered variable names
    center_vars   = ['pres', 'divv', 'ibmf', 'delp']
    face_vars     = ['velc', 'hvar']
    ins_vars      = ['velc', 'hvar', 'divv', 'pres', 'delp']
    poisson_vars  = ['delp', 'divv']
    imbound_vars  = ['ibmf', 'velc']

    scalar_info   = dict(tmax = 8, dt = 0.000625, Re = 100.0)

    simulation_info = dict(time_stepping = 'ab2', 
                           poisson_solver = 'superlu', 
                           maxiter = 2000,
                           tol = 1e-10,
                           with_ib = True,
                           mapping_type = 'classical')

    particle_info = [dict(input='HDF5', file='sm_body.00001.h5', vel = [0.0, -1.0])]

    # Define boundary conditions for variable pressure and velocity [left, right, bottom, top]
    bc_type_center = dict(delp = ['neumann', 'neumann', 'neumann', 'neumann'])
    bc_val_center  = dict(delp = [0.0, 0.0, 0.0, 0.0])

    bc_type_facex = dict(velc = ['dirichlet', 'dirichlet', 'dirichlet', 'dirichlet'])
    bc_val_facex  = dict(velc = [0.0, 0.0, 0.0, 0.0])

    bc_type_facey = dict(velc = ['dirichlet', 'dirichlet', 'dirichlet', 'dirichlet'])
    bc_val_facey  = dict(velc =[0.0, 0.0, 0.0, 0.0])

    # Create the grid and data
    gridc, gridx, gridy, scalars, particles = flowx.domain.Domain(nx, ny, xmin, xmax, ymin, ymax, \
                                                   center_vars, face_vars, scalar_info, particle_info, \
                                           bc_type_center=bc_type_center, bc_val_center=bc_val_center, \
                                               bc_type_facex=bc_type_facex, bc_val_facex=bc_val_facex, \
                                               bc_type_facey=bc_type_facey, bc_val_facey=bc_val_facey)

    domain_data_struct = [gridc, gridx, gridy, scalars, particles]

    poisson = flowx.poisson.Poisson(gridc, poisson_vars, poisson_info=simulation_info)

    imbound = flowx.imbound.ImBound(domain_data_struct, imbound_vars, imbound_info=simulation_info)

    ins = flowx.ins.IncompNS(poisson, imbound, domain_data_struct, ins_vars, ins_info=simulation_info)

    while scalars.time <= scalars.tmax:
    
        imbound.map_to_grid()

        ins.advance()

        for particle in particles: particle.advance()

        scalars.advance()
    
        # Display stats
        if scalars.nstep % 1 == 0: flowx.io.display_stats(scalars)  
        
if __name__ == "__main__":
    main()
