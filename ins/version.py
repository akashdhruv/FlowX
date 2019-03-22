"""Set up the version."""

import os


_version_major = 0
_version_minor = 1
_version_micro = ''
_version_extra = 'dev'

# construct full version string
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)
__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ['Development Status :: 1 - Alpha',
               'Environment :: Console',
               'License :: OSI Approved :: BSD 3-Clause License',
               'Operating System :: Unix',
               'Programming Language :: Python']

description = 'mae6225: incompressible Navier-Stokes solver'
# Long description will go up on the pypi page
long_description = """
mae6225
=====
mae6225 is a library that contains building blocks to solve the
two-dimensional incompressible Navier-Stokes equations using a fractional-step
method.

License
=======
``mae6225`` is licensed under the terms of the BSD 3-Clause license. See the
file "LICENSE" for information on the history of this software,
terms & conditions for usage, and a DISCLAIMER OF ALL WARRANTIES.
All trademarks referenced herein are property of their respective holders.
Copyright (c) 2019, Students of MAE-6225 (Spring 2019).
"""

NAME = 'mae6225'
MAINTAINER = ''
MAINTAINER_EMAIL = ''
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = 'https://github.com/Balaras-Group/MAE-6225'
DOWNLOAD_URL = ''
LICENSE = 'BSD 3-Clause'
AUTHOR = ''
AUTHOR_EMAIL = ''
PLATFORMS = 'Unix'
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGES = ['mae6225']
PACKAGE_DATA = {'mae6225': [os.path.join('styles', '*')]}
REQUIRES = ['numpy', 'matplotlib']
