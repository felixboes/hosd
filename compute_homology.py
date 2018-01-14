#!/usr/bin/env python

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
        raise "Python >= 2.7 is required for argument parsing."
    
    # Use the argparse library. The library optparse is deprecated since version 2.7.
    # Compare the documentation: https://docs.python.org/2/library/argparse.html
    
    # Create the argument parser.
    # Note: Config files can be processed. In order to do so, we have to give fromfile_prefix_chars='_' with _ a symbol.
    parser = argparse.ArgumentParser(
        add_help = True, 
        fromfile_prefix_chars='@', 
        description='Compute the homology of the compactification of the unilevel radial slit domains aka sulivan diagrams.'
    )
    
    # Provide all arguments.
    # Note: we provide nargs=N and the N arguments from the command line will be gathered together into a list.
    # Thus, we supress nargs.
    parser.add_argument('-g', '--gen',     required=True,  action='store', type=int, dest='g',              metavar='arg',  help='The genus of the Riemann surfaces')
    parser.add_argument('-m', '--pun',     required=True,  action='store', type=int, dest='m',              metavar='arg',  help='The number of punctures of the Riemann surfaces')
    parser.add_argument('-p', '--parametrized',            action='store_true',      dest='param_bdry',                     help='Consider parametrized boundary.', default=False)
    parser.add_argument('-v',                              action='store_true',      dest='more_verbose',                   help='Print more status information.', default=False)
    parser.add_argument('-s',                              action='store_true',      dest='silent',                         help='Print no status information.', default=False)
    parser.add_argument('--save_homchain', required=True,  action='store', nargs=1,  dest='homchain_file',  metavar='file', help='Store homchain file.')
    parser.add_argument('--only_good_stuff',               action='store_true',      dest='only_good_stuff',                help='Compute the homology of the quotient by the bad stuff.')
    args=vars( parser.parse_args() )

    # Setup verbosenes.
    if args['silent']:
        args['verbose'] = args.get('verbose', 'False')
        args['more_verbose'] = False
    elif args['more_verbose']:
        args['verbose'] = args.get('verbose', 'True')
    else:
        args['verbose'] = args.get('verbose', 'True')
    del args['silent']

    #Setup homchain
    args['homchain_file'] = args['homchain_file'][0]

    # Create directories.
    pyradbar.create_directories(args['g'], args['m'])
    
    # Tee the output to the result file.
    tee = pyradbar.Tee('./results/' + ''.join( [str(param).replace(' ', '_').replace('/', '_') for param in sys.argv if str(param) ] ).strip('.').lstrip('_'), 'w')

    # Write Preamble.
    pre, valid = pyradbar.preamble()
    sys.stdout.write(pre + '\n')
    sys.stdout.flush()
    
    # Stop of something went wrong.
    if valid == False:
        print "Could not initialize everything. Abroating."
        return 1
    
    # Start actual computation.
    pyradbar.compute_homology(**args )
    
    # Cleanup screen with new lines. This is not written to the result file.
    sys.stdout.write('\n\n\n')
    sys.stdout.flush()

if __name__ == "__main__":
    main()
