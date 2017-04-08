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
import subprocess
import sys

def call_sage(g=1, m=2, parametrized=False, ring=None, homol_alg=None, verbose=None, more_verbose=None, sanity_checks=None, only_good_stuff=None, sage_path=None, result_file=None, homcain_file=None ):
    script_path = './pyradbar/homology_computation.py'
    sys.stdout.write('Calling ' + sage_path + ' -python ' + script_path + ' ' + str(g) + ' ' + str(m) + ' ' + str(parametrized) + ' ' + str(ring) + ' ' + str(homol_alg) + ' ' + str(verbose) + ' ' + str(more_verbose) + ' ' + str(sanity_checks) + ' ' + str(only_good_stuff) + ' ' + str(result_file) + (' ' + str(homcain_file[0]) if homcain_file is not None else str('') ) + '\n')
    sys.stdout.flush()
    cmd = [sage_path, "-python", script_path, str(g), str(m), str(parametrized), str(ring), str(homol_alg), str(verbose), str(more_verbose), str(sanity_checks), str(only_good_stuff), str(result_file)]
    if homcain_file is not None:
        cmd.append(homcain_file[0])
    subprocess.call(cmd)

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
    parser.add_argument('-g', '--gen',   required=True,  action='store', type=int, dest='g',             metavar='arg',  help='The genus of the Riemann surfaces')
    parser.add_argument('-m', '--pun',   required=True,  action='store', type=int, dest='m',             metavar='arg',  help='The number of punctures of the Riemann surfaces')
    parser.add_argument('-r', '--coeff',                 action='store', type=str, dest='ring',          metavar='arg',  help='The coefficient ring. Either ZZ or a field in sage notation like QQ or GF(2).', default='QQ')
    parser.add_argument('-p', '--parametrized',          action='store_true',      dest='parametrized',                  help='Consider parametrized boundary.', default=False)
    parser.add_argument('--alg', '--homology_algorithm', action='store', type=str, dest='homol_alg',     metavar='arg',  help='The algorithm that is used to compute the homology e.g. \'chomp\'. All possible options are found in the sage manual.', default='auto')
    parser.add_argument('--sanity_checks',               action='store_true',      dest='sanity_checks',                 help='Perform sanity checks, e.g. D(n) \circ D(n+1) = 0.', default=False)
    parser.add_argument('-v',                            action='store_true',      dest='more_verbose',                  help='Print more status information.', default=False)
    parser.add_argument('-s',                            action='store_true',      dest='silent',                        help='Print no status information.', default=False)
    parser.add_argument('--sage',                        action='store', type=str, dest='sage_path',     metavar='path', help='The Path to the sage executable', default='./sage-6.8-x86_64-Linux/sage')
    parser.add_argument('--save_homchain',               action='store', nargs=1,  dest='homcain_file',  metavar='file', help='Store homchain file.')
    parser.add_argument('--only_good_stuff',             action='store_true',      dest='only_good_stuff',               help='Compute the homology of the quotient by the bad stuff.')
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

    # The name of the results file.
    args['result_file'] = './results/' + ''.join( [str(param).replace(' ', '_').replace('/', '_') for param in sys.argv if str(param) ] ).strip('.').lstrip('_')
    
    # Create directories.
    pyradbar.create_directories(args['g'], args['m'])
    
    # Tee the output to the result file.
    tee = pyradbar.Tee(args['result_file'], 'w')

    # Write Preamble.
    pre, valid = pyradbar.preamble( args['sage_path'] )
    sys.stdout.write(pre + '\n')
    sys.stdout.flush()
    
    # Since Sage uses file descriptors, we cannot use the Tee class directly.
    # That's unfortunate...
    tee.stop_logging()
    
    # Stop of something went wrong.
    if valid == False:
        print "Could not initialize everything. Abroating." 
        return 1
    
    # Start actual computation.
    call_sage(**args )
    
    # Cleanup screen with new lines. This is not written to the result file.
    sys.stdout.write('\n\n\n')
    sys.stdout.flush()

if __name__ == "__main__":
    main()
