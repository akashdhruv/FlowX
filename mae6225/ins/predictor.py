import numpy
from .convective_u import *
from .convective_v import *
from .diffusion_u  import *
from .diffusion_v  import *

def predictor(gridx, gridy, ivar, cvar, dvar, Re, dt):


    convective_u(gridx, gridy, ivar, cvar)
    convective_v(gridx, gridy, ivar, cvar)

    diffusion_u(gridx, ivar, dvar, Re)
    diffusion_v(gridy, ivar, dvar, Re)

    u  = gridx.get_values(ivar)
    v  = gridy.get_values(ivar)
    
    convx = gridx.get_values(cvar) 
    convy = gridy.get_values(cvar)

    diffx = gridx.get_values(dvar)
    diffy = gridy.get_values(dvar)

    u[1:-1,1:-1] = u[1:-1,1:-1] + dt*(convx[1:-1,1:-1] + diffx[1:-1,1:-1])
    v[1:-1,1:-1] = v[1:-1,1:-1] + dt*(convy[1:-1,1:-1] + diffy[1:-1,1:-1])

    gridx.fill_guard_cells(ivar)
    gridy.fill_guard_cells(ivar)

    return
