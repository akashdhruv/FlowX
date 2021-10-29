import numpy
import h5py

def write_hdf5_output(domain_data_struct=[None]*5, center_vars=None, face_vars=None):
    """
    Subroutine to write output to hdf5 file.

    Arguments
    ---------

    domain_data_struct : object list
    [gridc, gridx, gridy, scalars, particles]

    center_vars : string list
           List containing center var variable names

    face_vars : string list
          List containing face var variable names
    """

    gridc, gridx, gridy, scalars, particles = domain_data_struct

    if center_vars:
        for centervar in center_vars:
            map_center_to_edge()

    if face_vars:
        for facevar in face_vars:
            map_facex_to_edge()
            map_facey_to_edge()

    return


def map_center_to_edge(grid, ivar):

    return

def map_facex_to_edge(grid, ivar):

    return

def map_facey_to_edge(grid, ivar):

    return
