# Incompressible Navier-Stokes Solver

`mae6225` is a Python package that contains all necessary objects to build a solver for the two-dimensional incompressible Navier-Stokes equations.

`mae6225` is the result of a collaborative project between the students of the CFD class (MAE-6225) taught at the George Washington University (Spring 2019) by Prof. Balaras.

## Installation

---

### Dependencies (last tested)

* Python 3 (3.6.8)
* NumPy (1.15.4)
* Matplotlib (3.0.2)

To install the dependencies, we recommend using [`conda`](https://www.anaconda.com/distribution/).

To create a new `conda` environment to use `mae6225`:

```bash
conda create --name py36-mae6225 python=3.6 numpy matplotlib
conda activate py36-mae6225
```

Alternatively, you can use the file `environment.yaml`:

```bash
conda env create --file environment.yaml
conda activate py36-mae6225
```

### Installing `mae6225`

To install the Python package `mae6225`:

```bash
git clone https://github.com/Balaras-Group/MAE-6225.git
cd MAE-6225
python setup.py develop
```

### Run the tests

```bash
python tests/all.py
```

