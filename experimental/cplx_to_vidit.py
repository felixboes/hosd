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
import os
import time
import inspect

def call_sage(vidit_file_prefix, g=1, m=2, ring=None, verbose=None, sage_path=None):
    script_path = './pyradbar/cplx_to_vidit.py'
    sys.stdout.write('Calling ' + sage_path + ' -python ' + script_path + ' ' + vidit_file_prefix + ' ' + str(g) + ' ' + str(m) + ' ' + str(ring) + ' ' + str(verbose) + '\n')
    sys.stdout.flush()
    cmd = [sage_path, "-python", script_path, vidit_file_prefix, str(g), str(m), str(ring), str(verbose) ]
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
    parser.add_argument('--vidit',       required=True,  action='store', type=str,  dest='vidit_file_prefix',  metavar='file', help='Prefix for the file.')
    parser.add_argument('-g', '--gen',   required=True,  action='store', type=int, dest='g',             metavar='arg',  help='The genus of the Riemann surfaces')
    parser.add_argument('-m', '--pun',   required=True,  action='store', type=int, dest='m',             metavar='arg',  help='The number of punctures of the Riemann surfaces')
    parser.add_argument('-r', '--coeff',                 action='store', type=str, dest='ring',          metavar='arg',  help='The coefficient ring. Either ZZ or a field in sage notation like QQ or GF(2).', default='QQ')
    parser.add_argument('-s',                            action='store_true',      dest='silent',                        help='Print no status information.', default=False)
    parser.add_argument('--sage',                        action='store', type=str, dest='sage_path',     metavar='path', help='The Path to the sage executable', default='./sage-6.8-x86_64-Linux/sage')
    args=vars( parser.parse_args() )
    
    if args['silent']:
        args['verbose'] = args.get('verbose', 'False')
    else:
        args['verbose'] = args.get('verbose', 'True')
    del args['silent']

    pre, valid = pyradbar.preamble( args['sage_path'] )
    sys.stdout.write(pre + '\n')
    sys.stdout.flush()
    
    if valid == False:
        print "Could not initialize everything. Abroating." 
        return 1
    
    call_sage( **args )
    
    sys.stdout.write('\n\n\n')
    sys.stdout.flush()

if __name__ == "__main__":
    main()
