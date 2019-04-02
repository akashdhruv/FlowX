"""Routine to gather statistics."""

import numpy


def stats(gridc, gridx, gridy, ivar, pvar, divc):
    """Velocity correction in x and y direction.

    Arguments
    ---------
    gridc : grid object (cell center)
        Grid contaning data in cell center
    gridx : grid object (x-direction)
        Grid containing data in x-direction
    gridy : grid object  (y-direction)
        Grid containing data in y-direction
    ivar : string
        Name of the grid variable of the velocity solution
    pvar : string
        Name of the grid variable of the pressure solution
    divc : string
        Name of the grid variable for divergence

    Returns
    -------
    ins_stats : dictionary for stats

    """
    u = gridx.get_values(ivar)
    v = gridy.get_values(ivar)

    p = gridc.get_values(pvar)
    div = gridc.get_values(divc)

    max_u, min_u = numpy.amax(u), numpy.amin(u)
    max_v, min_v = numpy.amax(v), numpy.amin(v)
    max_p, min_p = numpy.amax(p), numpy.amin(p)
    max_div, min_div = numpy.amax(div), numpy.amin(div)

    ins_stats = dict()
    ins_stats['max_u'] = max_u
    ins_stats['min_u'] = min_u
    ins_stats['max_v'] = max_v
    ins_stats['min_v'] = min_v
    ins_stats['max_p'] = max_p
    ins_stats['min_p'] = min_p
    ins_stats['max_div'] = max_div
    ins_stats['min_div'] = min_div

    return ins_stats
