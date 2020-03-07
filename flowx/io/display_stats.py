"""Module with I/O functions."""

import numpy


def display_stats(scalars):
    """Display stats specified by user.

    Arguments
    ---------
    scalars : object
            Scalar class object containing neccesary information for head up display
    """
    print('------------ Time = {} ---------------'.format(scalars.time))
    print('Number of poisson iterations    : {}'.format(scalars.stats['ites']))
    print('Final poisson residual : {}'.format(scalars.stats['res']))
    print('Max, Min, U   : {}, {}'.format(scalars.stats['max_u'], scalars.stats['min_u']))
    print('Max, Min, V   : {}, {}'.format(scalars.stats['max_v'], scalars.stats['min_v']))
    print('Max, Min, P   : {}, {}'.format(scalars.stats['max_p'], scalars.stats['min_p']))
    print('Max, Min, DIV : {}, {}'.format(scalars.stats['max_div'], scalars.stats['min_div']))
    print('\n')
