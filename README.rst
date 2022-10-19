.. |icon| image:: ./media/icon.svg
  :width: 25 
 
============  
|icon| FlowX
============

|Code style: black|

|INS| |Poisson| |Grid| |Publish|

FlowX is a python based library to integrate Computational Fluid Dynamics (CFD) and geometry problems with quantum computing, and graphical libraries. It serves as sandbox to play with new features/algorithms and explore new research areas.

Installation
============

We recommend to install FlowX in development mode since it intended for development and design of new ideas. This can be easily accomplished using the ``setup`` script located in the project root directory and executing,

::

   ./setup develop

Development mode enables testing of features/updates directly from the source code and is an effective method for debugging. Note that the ``setup`` script relies on ``click``, which can be installed using,

::

  pip install click

Note that ``pip`` should point to ``python3+`` installation package ``pip3``, and FlowX has been tested on `python3.8+`

Testing
=======

Tests for basic functionality are located in ``tests/`` folder under project root directory and can be run using

.. code:: bash

   python3 tests/all.py
   
Examples
========

More interesting simulations involving fluid dynamics and quantum computing are located in ``examples/`` folder under project root directory, and is an evolving aspect of this project.

.. |Code style: black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   
.. |INS| image:: https://github.com/akashdhruv/FlowX/workflows/INS/badge.svg
.. |Poisson| image:: https://github.com/akashdhruv/FlowX/workflows/Poisson/badge.svg
.. |Grid| image:: https://github.com/akashdhruv/FlowX/workflows/Grid/badge.svg
.. |Publish| image:: https://github.com/akashdhruv/FlowX/workflows/Publish/badge.svg
