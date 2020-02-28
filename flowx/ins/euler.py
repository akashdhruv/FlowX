from .corrector import corrector
from .divergence import divergence
from .predictor import predictor
from .stats import stats

def euler(gridc, gridx, gridy, scalars, predcorr):

    if(predcorr == 'predictor'):

        # Calculate predicted velocity: u* = dt*H(u^n)
        predictor(gridx, gridy, 'velc', 'hvar', scalars.var['Re'], scalars.var['dt'])
    
        # Calculate RHS for the pressure Poission solver div(u)/dt
        divergence(gridc, gridx, gridy, 'velc', 'divp', ifac = scalars.var['dt'])


    elif(predcorr == 'corrector'):

        # Calculate corrected velocity u^n+1 = u* - dt * grad(P) 
        corrector(gridc, gridx, gridy, 'velc', 'pres', scalars.var['dt'])
    
        # Calculate divergence of the corrected velocity to display stats
        divergence(gridc, gridx, gridy, 'velc', 'divc')
    
        # Calculate stats
        scalars.stats.update(stats(gridc, gridx, gridy, 'velc', 'pres', 'divc'))
