"""User defined module for simulation."""
  
import numpy
import flowx.archive
import simulation

def example():
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

    simulation.update_bc_val(gridx, gridy, 'velc', scalars.to)
    simulation.get_initial(gridc, gridx, gridy, 'velc', 'pres')

    while scalars.time <= scalars.tmax:
    
        # Update the time-dependent boundary condition value
        simulation.update_bc_val(gridx, gridy, 'velc', scalars.time)
    
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
    example()
