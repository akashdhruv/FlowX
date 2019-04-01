"""Module with I/O functions."""

import numpy

def display_stats(t,ites,res,stats):
    """display stats specified by user

    Arguments
    ---------

    """

    print('-----Time = {}-------'.format(t))
    print('Number of poisson iterations: {}'.format(ites))
    print('Final poisson residual residual: {}'.format(res))
    print('Max, Min, U: {} {}'.format(stats[0],stats[1]))
    print('Max, Min, V: {} {}'.format(stats[2],stats[3]))
    print('Max, Min, P: {} {}'.format(stats[4],stats[5]))
    print('Max, Min, DIV: {} {}'.format(stats[6],stats[7]))
    print('\n')
