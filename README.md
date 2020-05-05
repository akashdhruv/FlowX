# Incompressible Navier-Stokes Solver

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/Balaras-Group/flowX/raw/master/LICENSE)
[![Travis](https://img.shields.io/travis/com/Balaras-Group/flowX/develop.svg?style=flat-square&logo=travis)](https://travis-ci.com/Balaras-Group/flowX)

`flowx` is a Python package that contains all necessary objects to build a solver for the two-dimensional incompressible Navier-Stokes equations.

`flowx` is the result of a collaborative project between the students of the CFD class (MAE-6225) taught at the George Washington University.

## Installation

---

### Dependencies (last tested)

* Python 3 (3.6.8)
* NumPy (1.15.4)
* Matplotlib (3.0.2)

To install the dependencies, we recommend using [`conda`](https://www.anaconda.com/distribution/).

To create a new `conda` environment to use `flowx`:

```bash
conda create --name py36-flowx python=3.6 numpy matplotlib
conda activate py36-flowx
```

Alternatively, you can use the file `environment.yaml`:

```bash
conda env create --file environment.yaml
conda activate py36-flowx
```

### Installing `flowx`

To install the Python package `flowx`:

```bash
git clone https://github.com/Balaras-Group/flowX.git
cd flowX
python setup.py develop
```

### Run the tests

```bash
python tests/all.py
```
