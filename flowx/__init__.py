import argparse,os

parser = argparse.ArgumentParser()
parser.add_argument("--env", help="argument for serial/parallel environment", default='serial')
parser.add_argument("-f", "--fff", help="a dummy argument to fool ipython", default="1")

args = parser.parse_args().__dict__

__environment__ = args['env']

if __environment__ in ['serial','SERIAL']:
    __environment__ = 'serial'

elif __environment__ in ['parallel','PARALLEL']:
    __environment__ = 'parallel'

else:
    raise ValueError("".join(['Environment:',__environment__,' not recognized']))

from flowx.version import __version__
from flowx.domain import domain_main
from flowx.poisson import poisson_main
from flowx.ins import ins_main
from flowx.imbound import imbound_main
from flowx.quantum import quantum_main
import flowx.io
