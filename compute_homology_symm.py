#!/usr/bin/env python3

# The software pyradbar is a bunch of programs to compute the homology of
# Sullivan diagrams.
# Copyright (C) 2015 - 2017  Felix Boes
#
# This file is part of pyradbar.
#
# pyradbar is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyradbar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyradbar.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import pyradbar
import sys

def main():
    # check for correct version.
    major, minor = sys.version_info[0], sys.version_info[1]
    if major < 2 or (major == 2 and minor < 7):
        raise RuntimeError('Python >= 2.7 is required for argument parsing.')

    # Use the argparse library. The library optparse is deprecated since version 2.7.
    # Compare the documentation: https://docs.python.org/2/library/argparse.html

    # Create the argument parser.
    # Note: Config files can be processed. In order to do so, we have to give fromfile_prefix_chars='_' with _ a symbol.
    parser = argparse.ArgumentParser(
        add_help = True,
        fromfile_prefix_chars='@',
        description='Compute the homology of Symm_d in a range.'
    )

    # Provide all arguments.
    # Note: we provide nargs=N and the N arguments from the command line will be gathered together into a list.
    # Thus, we supress nargs.
    parser.add_argument('-d', '--max_d', required=True,  action='store', type=int, dest='max_d',         metavar='arg',  help='The top degree of the Symm_d complex.')
    parser.add_argument('-c', '--num_cyc',               action='store', type=int, dest='num_cyc',       metavar='arg',  help='Restrict to a specific number of cycles.', default=0)
    parser.add_argument('-v',                            action='store_true',      dest='verbose',                       help='Print status information.', default=False)
    parser.add_argument('-s',                            action='store_true',      dest='silent',                        help='Print no status information.', default=False)
    parser.add_argument('--save_homchain',               action='store', nargs=1,  dest='homchain_file', metavar='file', help='Store homchain file.')
    args=vars( parser.parse_args() )

    # Setup verbosenes.
    if args['num_cyc'] <= 0:
      del args['num_cyc']
    args['homchain_file'] = args['homchain_file'][0]
    if args['silent']:
        args['verbose'] = args.get('verbose', 'False')
    else:
        args['verbose'] = args.get('verbose', 'True')
    del args['silent']


    # The name of the results file.
    result_file = './results/' + ''.join( [str(param).replace(' ', '_').replace('/', '_') for param in sys.argv if str(param) ] ).strip('.').lstrip('_')

    # Tee the output to the result file.
    tee = pyradbar.Tee(result_file, 'w')

    # Write Preamble.
    pre, valid = pyradbar.preamble()
    print(pre)

    # Stop of something went wrong.
    if valid == False:
        sys.stdout.write('Could not initialize everything. Abroating.\n')
        return 1

    # Start actual computation.
    pyradbar.compute_homology_symm(**args)

    # Cleanup screen with new lines.
    sys.stdout.write('\n')

if __name__ == '__main__':
    main()
