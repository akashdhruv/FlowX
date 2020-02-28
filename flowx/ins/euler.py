from .projection import predictor, corrector, divergence
from .stats import stats

def euler(gridc, gridx, gridy, scalars, var_list, predcorr):

    velc = var_list[0]
    hvar = var_list[1]
    divv = var_list[2]
    pres = var_list[3]

    if(predcorr == 'predictor'):

        # Calculate predicted velocity: u* = dt*H(u^n)
        predictor(gridx, gridy, velc, hvar, scalars.var['Re'], scalars.var['dt'])
    
        # Calculate RHS for the pressure Poission solver div(u)/dt
        divergence(gridc, gridx, gridy, velc, divv, ifac = scalars.var['dt'])


    elif(predcorr == 'corrector'):

        # Calculate corrected velocity u^n+1 = u* - dt * grad(P) 
        corrector(gridc, gridx, gridy, velc, pres, scalars.var['dt'])
    
        # Calculate divergence of the corrected velocity to display stats
        divergence(gridc, gridx, gridy, velc, divv)
    
        # Calculate stats
        scalars.stats.update(stats(gridc, gridx, gridy, velc, pres, divv))
