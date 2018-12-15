# generate_hte10_lattice

Constructs an input file for the high-temperature series expansion code [HTE10](http://wasd.urz.uni-magdeburg.de/jschulen/HTE10/) by A. Lohmann and J. Richter, see [Phys. Rev. B **89**, 014415 (2014)](https://doi.org/10.1103/PhysRevB.89.014415).

Requirements
------------

The code requires a python interpreter (2.7 or >=3.5) and NumPy (>=1.10).

How to run
----------

1. Create an input file describing the magnetic exchanges in the unit cell of a periodic spin lattice. Exchanges couple the *i*th spin in the unit cell with the *j*th spin, which can be located either in the same cell or one of the neighboring cells. This file has the following structure:

* Column 1: index of the exchange (**starting from 0**);
* Column 2: *i*th spin in the (0 0 0) unit cell;
* Column 3: *j*th spin;
* Columns 4, 5, 6: the index of the cell accommodating the *j*th spin, assuming that the *i*th spin lives in the (0 0 0) cell.

Example: for the kagome lattice, all exchanges are equivalent (zeros in the first column) we have three spins in the unit cell (indices 0, 1, and 2 in the second and thirds columns), and the input file `kagome.def` is:

```
0 0 1  0  0  0
0 0 2  0  0  0
0 1 2  0  0  0
0 0 1 -1  0  0
0 1 2  1 -1  0
0 0 2  0 -1  0
```

2. Decide how many cells to print along the three spatial dimensions, keeping in mind that the **boundary conditions are always periodic.** If a model is 2D in the *xy* plane, the last index should be unity; for 1D models (along *x*), two last indices should be unity.

3. Execute the script. For instance, for a kagome lattice with 20 by 20 cells:

```
./generate_hte10_lattice.py kagome.def -c 20 20 1
```

Usage
-----

```
./generate_hte10_lattice.py --help
```
