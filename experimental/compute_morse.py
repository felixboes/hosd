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

def call_sage(g=1, m=2, more_verbose=None, result_file=None, sage_path=None):
    script_path = './pyradbar/morse_computation.py'
    sys.stdout.write('Calling ' + sage_path + ' -python ' + script_path + ' ' + str(g) + ' ' + str(m) + ' ' + str(more_verbose) + ' ' + str(result_file) + '\n')
    sys.stdout.flush()
    cmd = [sage_path, "-python", script_path, str(g), str(m), str(more_verbose), str(result_file)]
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
    parser.add_argument('-v',                            action='store_true',      dest='more_verbose',                  help='Print more status information.', default=False)
    parser.add_argument('--sage',                        action='store', type=str, dest='sage_path',     metavar='path', help='The Path to the sage executable', default='./sage-6.8-x86_64-Linux/sage')
    args=vars( parser.parse_args() )

    # The name of the results file.
    args['result_file'] = './results/' + ''.join( [str(param).replace(' ', '_') for param in sys.argv if str(param) ] )

    tee = pyradbar.Tee(args['result_file'], 'w')

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
