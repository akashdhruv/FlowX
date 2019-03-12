import numpy

def convective_u(gridx, gridy, ivar, cvar):
    idxX = gridx.get_variable_indices(ivar)
    idxY = gridy.get_variable_indices(ivar)
    idcX = gridx.get_variable_indices(cvar)

    veluE = gridx.data[1:-1, 2:, idxX]
    veluW = gridx.data[1:-1,:-2, idxX]
    
    velvN = gridy.data[2:, 1:-1, idxY]
    velvS = gridy.data[:-2,1:-1, idxY]
    velvE = gridy.data[1:-1, 2:, idxY]
    velvW = gridy.data[1:-1,:-2, idxY]

    veluP = (veluE + veluW) / 2
    
    velvNW = (velvN + velvW) / 2
    velvNE = (velvN + velvE) / 2
    velvSW = (velvS + velvW) / 2
    velvSE = (velvS + velvE) / 2

    F[1:-1, 1:-1] = (-((veluP + veluE)**2 - (veluP veluW)**2) / 4 * gridx.dx - 
                     ((veluP + veluN)*(velvNE + velvNW)) / 4 * gridy.dy - 
                     ((veluP + veluS)*(velvSE + velvSW)) / 4 * gridy.dy)
    
    gridx.data[1:-1, 1:-1, idcX] = F[1:-1, 1:-1]
    
    return
