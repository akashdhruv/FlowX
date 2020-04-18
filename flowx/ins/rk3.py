"""Runge and Kutta 3 stage explicit time advancement routine"""

from .projection import predictor_rk3, corrector_rk3, divergence_rk3
from .stats import stats

def rk3(gridc, gridx, gridy, scalars, grid_var_list, predcorr, i=0):

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

    predcorr : string
           Flag for the fractional step method equations - 'predictor', 'divergence', 'corrector'

    i : integer 
           Number of stage

    """

    velc = grid_var_list[0]
    hvar = grid_var_list[1]
    divv = grid_var_list[2]
    pres = grid_var_list[3]

    if(predcorr == 'predictor'):

        # Calculate predicted velocity: u* = dt*H(u^n)
        predictor_rk3(gridx, gridy, velc, hvar, scalars.variable['Re'], scalars.variable['dt'], i=0)


    if(predcorr == 'divergence'):    
        # Calculate RHS for the pressure Poission solver div(u)/dt
        divergence_rk3(gridc, gridx, gridy, velc, divv, ifac = scalars.variable['dt'], i=0)


    elif(predcorr == 'corrector'):

        # Calculate corrected velocity u^n+1 = u* - dt * grad(P) 
        corrector_rk3(gridc, gridx, gridy, velc, pres, scalars.variable['dt'], i=0)
    
        # Calculate divergence of the corrected velocity to display stats
        divergence_rk3(gridc, gridx, gridy, velc, divv, i=0)
    
        # Calculate stats
        scalars.stats.update(stats(gridc, gridx, gridy, velc, pres, divv))
