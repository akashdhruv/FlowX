#Diffusion function
import numpy
def diffusion_v(gridy,vvar,dfvar, Re):
    
    #Define input variables
    #Re = Reynold's number
    #gridx = grid for x coordinates and values
    #gridy = grid for y coordinates and values
    #vij = y velocity at position i,j
    
    #call variables
    i_vvar= gridy.get_variable_indices(vvar)  #store position of velocity values
    nx, ny = gridy.nx, gridy.ny
    dx, dy = gridy.dx, gridy.dy  #store grid size in x and y dimensions 

    #calculate and store diffusion values
    nu = 1/Re # viscosity term as a function of Reynold's number
    
    # loop to populate duffion values at grid points
    D = numpy.zeros((nx+2,ny+1))
    
    for j in range(1,ny):
        for i in range(1,nx+1):
              D[i,j] = nu *(gridy.data[i+1,j,i_vvar] - 2*gridy.data[i,j,i_vvar] + gridy.data[i-1,j,i_vvar]) / dx**2 +                                    (gridy.data[i,j+1,i_ivar] - 2*gridy.data[i,j,i_vvar] + gridy.data[i,j-1, i_vvar]) / dy**2
    
    #Store diffusion values in grid object
    gridy.set_values(dfvar,D)
    return  # nothing to return since diffusion variable is stored
