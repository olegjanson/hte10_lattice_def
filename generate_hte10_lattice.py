#!/usr/bin/env python
"""Constructs an input file for the HTSE code by A. Lohmann and J. Richter."""

__author__ = "Oleg Janson"
__date__ = "Feb 20, 2020"
__email__ = "olegjanson@gmail.com"
__version__ = "1.0"


import argparse as ap
import textwrap
import numpy as np

DESCRIPTION = ("Constructs an input file for the HTSE code " +
               "by A. Lohmann and J. Richter.")
EPILOG = "Please report any bugs to Oleg Janson <olegjanson@gmail.com>."


def get_exchange_name(j_ij_index):
    '''Converts index of an exchange into the name: 0 --> "j1", 1 --> "j2" etc.

    Args:
        j_ij_index (int): Index of the exchange (starting from 0)

    Returns:
        (str): The name of the exchange.

    '''
    return "j{:d}".format(1 + j_ij_index)


def apply_pbc(base_vector, increment_vector, boundaries):
    '''Adds two vectors subject to periodic boundary conditions (PBC).

    Args:
        base_vector (np.array): The first vector.
        increment (np.array): The second vector.
        boundary (np.array): The maximal value of the sum (also a vector).

    Returns:
        (np.array): The sum of two vectors subject to PBC.

    '''
    return (base_vector + increment_vector) % boundaries


def get_index(spin, pos, cells):
    '''Returns the index of a certain spin in the spin lattice.

    '''
    return spin + nspins*(pos[0] + pos[1]*cells[0] + pos[2]*cells[0]*cells[1])


if __name__ == "__main__":
    parser = ap.ArgumentParser(description=DESCRIPTION,
                               conflict_handler="resolve",
                               formatter_class=ap.RawTextHelpFormatter,
                               epilog=EPILOG)
    parser.add_argument('filename', help=textwrap.dedent('''\
the name of a six-column file with integer entries.\

Column 1: index of the respective exchange (starting from 0);
Column 2: \"i\"th spin (in the original unit cell);
Column 3: \"j\"th spin;
Columns 4, 5, 6: the index of the cell accommodating the \"j\"th spin,
                 assuming that the \"i\"th spin lives in the (0 0 0) cell.

Example: for the kagome lattice, the file is

0 0 1  0  0  0
0 0 2  0  0  0
0 1 2  0  0  0
0 0 1 -1  0  0
0 1 2  1 -1  0
0 0 2  0 -1  0
    '''))

    parser.add_argument("-c", "--cells", nargs=3, type=int, required=True,
                        help="the number of cells along three dimensions")
    parser.add_argument("-l", "--lattice", default='periodic lattice',
                        help="lattice name (Default: \"periodic lattice\")")
    parser.add_argument("-v", "--version", action="version",
                        help="print the version",
                        version="%(prog)s version {:s}".format(__version__))
    args = parser.parse_args()

    filename, cells, lattice = args.filename, tuple(args.cells), args.lattice

    spin_model = np.loadtxt(filename, dtype=np.int32)
    js = [dict(zip(["j_ij_index", "s_i", "s_j"], j_ij))
          for j_ij in spin_model[:, :3]]
    for _, j in enumerate(js):
        j['r_ij'] = np.array(spin_model[_, 3:])

    ncells = np.prod(cells)
    nspins = 1 + np.max(spin_model[..., 1:3])  # max spin index

    out = (("# {:s} N={:d}, lattice {:d}x{:d}x{:d},"
            + "periodic boundary conditions\n")
           .format(lattice, ncells*nspins, *cells))
    out += ("# Number of sites | Number of bonds "
            + "| Number of sites in the unit cell\n")
    out += ("{:3d} {:3d} {:2d}\n"
            .format(ncells*nspins, ncells*spin_model.shape[0], nspins))
    out += "# the numbers of the sites in the central unit cell\n"
    out += "\n".join("{:3d}".format(spin) for spin in range(nspins))
    out += "\n# Bond s1 s2\n"

    htse = []
    all_js = np.arange(ncells * len(js)).reshape(cells + (len(js),))
    c = np.array(cells)
    for ((x, y, z, _j), _) in np.ndenumerate(all_js):
        r, j = np.array([x, y, z]), js[_j]
        htse.append(" {:2d} {:2d} {:2d} {:s}\n"
                    .format(_, get_index(j["s_i"], r, c),
                            get_index(j["s_j"], apply_pbc(j["r_ij"], r, c), c),
                            get_exchange_name(j["j_ij_index"])))

    out += "".join(htse)
    out += "# end of file"
    print(out)
