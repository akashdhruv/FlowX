import numpy
import flowx

def main():

    # Define grid parameters
    nx, ny = 40, 40
    xmin, xmax = -0.5, 0.5
    ymin, ymax = -0.5, 0.5

    # Define cell-centered variable names
    center_vars   = ['pres', 'divv', 'ibmf', 'ibmx', 'ibmy', 'delp']
    face_vars     = ['velc', 'hvar']
    ins_vars      = ['velc', 'hvar', 'divv', 'pres', 'delp']
    poisson_vars  = ['delp', 'divv']
    imbound_vars  = ['ibmf', 'velc', 'ibmx', 'ibmy']

    scalar_info   = dict(tmax = 20, dt = 0.001, Re = 100.0, Re_s = 10.0, mu_s = 1.0)

    simulation_info = dict(with_ib = True, ib_type = 'visco', extrap_solid = 10)

    particle_info = [dict(input='HDF5', file='sm_body.00001.h5')]

    # Define boundary conditions for variable pressure and velocity [left, right, bottom, top]
    bc_type_center = dict(delp = ['neumann', 'neumann', 'neumann', 'neumann'],
                          ibmf = ['projection', 'projection', 'projection', 'projection'],
                          ibmx = ['projection', 'projection', 'projection', 'projection'],
                          ibmy = ['projection', 'projection', 'projection', 'projection'])

    bc_val_center  = dict(delp = [0.0, 0.0, 0.0, 0.0])

    bc_type_facex = dict(velc = ['dirichlet', 'dirichlet', 'dirichlet', 'dirichlet'])
    bc_val_facex  = dict(velc = [0.0, 0.0, 0.0, 1.0])

    bc_type_facey = dict(velc = ['dirichlet', 'dirichlet', 'dirichlet', 'dirichlet'])
    bc_val_facey  = dict(velc = [0.0, 0.0, 0.0, 0.0])

    # Create the grid and data
    gridc, gridx, gridy, scalars, particles = flowx.serial.domain_main(nx, ny, xmin, xmax, ymin, ymax, 
                                              center_vars, face_vars, scalar_info, particle_info, 
                                              bc_type_center=bc_type_center, bc_val_center=bc_val_center, 
                                              bc_type_facex=bc_type_facex, bc_val_facex=bc_val_facex, 
                                              bc_type_facey=bc_type_facey, bc_val_facey=bc_val_facey)

    domain_data_struct = [gridc, gridx, gridy, scalars, particles]

    poisson = flowx.poisson_main(gridc, poisson_vars, poisson_info=simulation_info)

    imbound = flowx.imbound_main(domain_data_struct, imbound_vars, imbound_info=simulation_info)

    ins = flowx.ins_main(poisson, imbound, domain_data_struct, ins_vars, ins_info=simulation_info)

    imbound.map_to_grid()

    while scalars.time <= scalars.tmax:
    
        ins.advance()
    
        imbound.advect()
    
        # Display stats
        if scalars.nstep % 1 == 0:
            print("Level Set Advection Time: ",imbound._advection_time)
            flowx.io.display_stats(scalars) 
        
        #if scalars.nstep % 100 == 0:
        #    flowx.io.plot_contour_zero(gridc, scalars, 'ibmf', 'ibmx', 'ibmy')


        scalars.advance()

if __name__ == "__main__":
    main()
