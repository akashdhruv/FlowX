"""Module with I/O functions."""

import numpy


def display_stats(t, stats):
    """Display stats specified by user.

    Arguments
    ---------
    t : float
        The time.
    stats : dictionary
        Dictionary with the stats to display.

    """
    print('------------ Time = {} ---------------'.format(t))
    print('Number of poisson iterations    : {}'.format(stats['ites']))
    print('Final poisson residual : {}'.format(stats['res']))
    print('Max, Min, U   : {}, {}'.format(stats['max_u'], stats['min_u']))
    print('Max, Min, V   : {}, {}'.format(stats['max_v'], stats['min_v']))
    print('Max, Min, P   : {}, {}'.format(stats['max_p'], stats['min_p']))
    print('Max, Min, DIV : {}, {}'.format(stats['max_div'], stats['min_div']))
    print('\n')
