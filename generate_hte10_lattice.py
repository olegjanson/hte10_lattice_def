#!/usr/bin/env python
""" Constructs an input file for the HTSE code by A. Lohmann and J. Richter. """

__author__  = 'Oleg Janson'
__date__    = 'Dec 15, 2018'
__email__   = 'olegjanson@gmail.com'
__version__ = '0.93'

import argparse
import textwrap
import numpy as np

description_text = \
    "Constructs an input file for the HTSE code by A. Lohmann and J. Richter."
epilog_text = "Please report any bugs to Oleg Janson <olegjanson@gmail.com>."

parser = argparse.ArgumentParser(description = description_text,
                                 conflict_handler = "resolve",
                                 formatter_class =\
                                     argparse.RawTextHelpFormatter,
                                 epilog = epilog_text)
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

parser.add_argument('-c', '--cells', nargs = 3, type = int, required = True,
                    help = 'the number of cells along three dimensions')
parser.add_argument('-l', '--lattice', default='some lattice',\
                    help='the name of the lattice (Default: \"some lattice\")')
parser.add_argument('-v', '--version', action='version',\
                     help='print the version',\
                     version='%(prog)s version {:s}'.format(__version__))
args = parser.parse_args()

filename, cells, lattice = args.filename, tuple(args.cells), args.lattice
cells_along_x, cells_along_y, cells_along_z = cells

spin_model_as_array = np.loadtxt(filename, dtype=np.int32)
exchanges = [dict(zip(['index','spin_i','spin_j','Tx','Ty','Tz'], j_ij)) \
             for j_ij in spin_model_as_array]
number_of_exchanges = spin_model_as_array.shape[0]
number_of_cells = np.prod(cells)
number_of_spins = 1 + np.max(spin_model_as_array[..., 1:3]) # max spin index

def get_exchange_name(exchange_index):
    ''' Returns the name of the exchange: 0 --> j1, 1 --> j2 etc. '''
    return 'j{:d}'.format(1 + exchange_index)

def apply_pbc(component_of_a_vector, increment, boundary):
    ''' Enforces 1D periodic boundary conditions. '''
    return (component_of_a_vector + increment)  %  boundary

def get_spin_index(spin, number_of_spins, x, y, z,
                   cells_along_x, cells_along_y):
    ''' Returns the index of a certain spin in the spin lattice. '''
    return spin  +  x * number_of_spins \
                 +  y * number_of_spins * cells_along_x \
                 +  z * number_of_spins * cells_along_x * cells_along_y

out = (  '# {:s} N={:d}, lattice {:d}x{:d}x{:d},' \
       + 'periodic boundary conditions\n')\
      .format(lattice, number_of_cells * number_of_spins, *cells)
out += '# Number of sites | Number of bonds ' \
    +  '| Number of sites in the unit cell\n'
out += '{:3d} {:3d} {:2d}\n'.format(number_of_cells * number_of_spins,
                                    number_of_cells * number_of_exchanges,
                                    number_of_spins)
out += '# the numbers of the sites in the central unit cell\n'
out += '\n'.join('{:3d}'.format(spin) for spin in range(number_of_spins))
out += '\n# Bond s1 s2\n'

htse_def = []
for x, y, z in np.ndindex((cells)):
    spin_index = x * cells_along_y * cells_along_z  +  y * cells_along_z  +  z
    spin_index *= number_of_exchanges
    for exchange_index, exchange in enumerate(exchanges):
        htse_def.append(' {:2d}'.format(spin_index + exchange_index))
        htse_def.append(' {:2d}'.format(get_spin_index(exchange['spin_i'],
                                                       number_of_spins,
                                                       x,
                                                       y,
                                                       z,
                                                       cells_along_x,
                                                       cells_along_y)))
        htse_def.append(' {:2d}'.format(get_spin_index(exchange['spin_j'],
                                                       number_of_spins,
                                                       apply_pbc(
                                                           exchange['Tx'],
                                                           x,
                                                           cells_along_x),
                                                       apply_pbc(
                                                           exchange['Ty'],
                                                           y,
                                                           cells_along_y),
                                                       apply_pbc(
                                                           exchange['Tz'],
                                                           z,
                                                           cells_along_z),
                                                       cells_along_x,
                                                       cells_along_y)))
        htse_def.append(' {:s}\n'.format(get_exchange_name(exchange['index'])))

out += ''.join(htse_def)
out += '# end of file'
print(out)
