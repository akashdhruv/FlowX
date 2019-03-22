import numpy

def convective_u(gridx, gridy, ivar, cvar):
    idxX = gridx.get_variable_indices(ivar)
    idxY = gridy.get_variable_indices(ivar)
    idcX = gridx.get_variable_indices(cvar)

    veluE = gridx.data[2:,  1:-1, idxX]
    veluW = gridx.data[:-2, 1:-1, idxX]
    veluP = gridx.data[1:-1,1:-1, idxX]
    veluN = gridx.data[1:-1,  2:, idxX]
    veluS = gridx.data[1:-1, :-2, idxX]
   
    velvNW = gridy.data[1:-2, 1:, idxY]
    velvNE = gridy.data[2:-1, 1:, idxY]
    velvSW = gridy.data[1:-2,:-1, idxY]
    velvSE = gridy.data[2:-1,:-1, idxY]
    
    F = (-((veluP + veluE)**2 - (veluP + veluW)**2) / (4 * gridx.dx) - 
          ((veluP + veluN)*(velvNE + velvNW)) / (4 * gridy.dy) - 
          ((veluP + veluS)*(velvSE + velvSW)) / (4 * gridy.dy))
    
    gridx.data[1:-1, 1:-1, idcX] = F
    
    return
