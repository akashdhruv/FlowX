"""Euler explicit time advancement routine"""

from flowx.poisson.poisson_main import solve_poisson
from flowx.ins.solvers.projection import predictor, predictor_AB2, predictor_RK3, corrector, divergence
from flowx.ins.solvers.stats import stats
from flowx.ins.solvers.mass_balance import get_qin, get_qout, rescale_velocity, get_convvel, update_outflow_bc
from flowx.imbound.imbound_main import solve_imbound

def ins_advance(gridc, gridx, gridy, scalars, grid_var_list, **kwargs):

    """
    Subroutine for the fractional step euler explicit time advancement of Navier Stokes equations
 
    Arguments
    ---------
    gridc : object
          Grid object for cell centered variables

    gridx : object
          Grid object for x-face variables

    gridy : object
          Grid object for y-face variables

    scalars: object
           Scalars object to access time-step and Reynold number

    grid_var_list : list
           List containing variable names for velocity, RHS term from the previous time-step, divergence and pressure

    **kwargs : list
             List of optional keywords arguments for time-stepping and poisson solver

             Accepted keywords and their options

             time-stepping = 'euler', 'ab2', 'rk3' -> default 'ab2'
    """

    _time_stepping = 'ab2'
    _with_ib = False

    if 'time_stepping' in kwargs: _time_stepping = kwargs.get('time_stepping')
    if 'with_ib' in kwargs: _with_ib = kwargs.get('with_ib')

    if _time_stepping is 'euler':
        solve_predictor = predictor
    elif _time_stepping is 'ab2':
        solve_predictor = predictor_AB2
    elif _time_stepping is 'rk3':
        solve_predictor = predictor_RK3

    velc = grid_var_list[0]
    hvar = grid_var_list[1]
    divv = grid_var_list[2]
    pres = grid_var_list[3]

    # Compute mass in
    Qin =  get_qin(gridx, velc) + get_qin(gridy, velc)

    # Update BC for predictor step
    update_outflow_bc(gridx, velc, scalars.variable['dt'])
    update_outflow_bc(gridy, velc, scalars.variable['dt'])

    # Calculate predicted velocity: u* = dt*H(u^n)       
    solve_predictor(gridx, gridy, velc, hvar, scalars.variable['Re'], scalars.variable['dt'])
 
    # Immersed boundary forcing
    if _with_ib: solve_imbound()

    gridx.fill_guard_cells(velc)
    gridy.fill_guard_cells(velc)

    # Calculate RHS for the pressure Poission solver div(u)/dt
    divergence(gridc, gridx, gridy, velc, divv, ifac = scalars.variable['dt'])
    gridc.fill_guard_cells(divv)

    # Compute mass out
    Qout =  get_qout(gridx, velc) + get_qout(gridy, velc)

    # Rescale velocity to balance mass
    rescale_velocity(gridx, velc, Qin, Qout)
    rescale_velocity(gridy, velc, Qin, Qout)

    # Update BC for corrector step
    update_outflow_bc(gridx, velc, scalars.variable['dt'], convvel=[0.0,0.0,0.0,0.0])
    update_outflow_bc(gridy, velc, scalars.variable['dt'], convvel=[0.0,0.0,0.0,0.0])

    # Solve pressure Poisson equation
    scalars.stats['ites'], scalars.stats['res'] = solve_poisson(gridc, pres, divv, **kwargs)

    # Calculate corrected velocity u^n+1 = u* - dt * grad(P) 
    corrector(gridc, gridx, gridy, velc, pres, scalars.variable['dt'])
    gridx.fill_guard_cells(velc)
    gridy.fill_guard_cells(velc)
   
    # Calculate divergence of the corrected velocity to display stats
    divergence(gridc, gridx, gridy, velc, divv)
    gridc.fill_guard_cells(divv)

    # Calculate stats
    scalars.stats.update(stats(gridc, gridx, gridy, velc, pres, divv))
