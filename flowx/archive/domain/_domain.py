from . import Grid
from . import Scalars
from . import Particles

def Domain(nx, ny, xmin, xmax, ymin, ymax, \
           center_vars=None, face_vars=None, scalar_info=None, particle_info=None, \
           bc_type_center=None, bc_val_center=None, \
           bc_type_facex=None, bc_val_facex=None, \
           bc_type_facey=None, bc_val_facey=None):


    gridc, gridx, gridy, scalars, particles = [object] * 5   
    
    if center_vars: 
        gridc = Grid('cell-centered', center_vars, nx, ny, xmin, xmax, ymin, ymax, \
                      user_bc_type=bc_type_center, user_bc_val=bc_val_center)

    if face_vars:
        gridx = Grid('x-face', face_vars, nx, ny, xmin, xmax, ymin, ymax, \
                      user_bc_type=bc_type_facex, user_bc_val=bc_val_facex)

        gridy = Grid('y-face', face_vars, nx, ny, xmin, xmax, ymin, ymax, \
                      user_bc_type=bc_type_facey, user_bc_val=bc_val_facey)

    if scalar_info: 
        scalars = Scalars(scalar_info)

    if particle_info:
        particles = [Particles(info, xmin, xmax, ymin, ymax, scalars) 
                     for info in particle_info]

    domain_data_struct = [gridc, gridx, gridy, scalars, particles] 
  
    return domain_data_struct
