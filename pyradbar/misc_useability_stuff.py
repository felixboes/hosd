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

import os
import platform
import subprocess
import sys
import time

def create_directories( g, m ):
    try:
        if not os.path.exists('./data/'):
            os.makedirs('./data/')
        if not os.path.exists('./results/'):
            os.makedirs('./results/')
    except:
        raise


def preamble( sage_path = None ):
    """"Print a nice preamble."""
    ret  = 'Programcall:     '
    for arg in sys.argv:
        ret += ' ' + arg
    ret += '\n'
    
    version_by_git = None
    found_git = False
    sage_version = None
    found_sage = False
    
    # Find and store git version if available.
    try:
        if os.path.exists('.git'):
            try:
                version_by_git = subprocess.check_output( [ 'git', 'rev-parse', 'HEAD' ], stderr=open(os.devnull, 'wb') ).replace('\n', '')
                found_git = True
            except:
                # Calling git returned an error. We ignore this.
                pass
    except:
        # Looking up the existence of a .git directory raised an exception.
        raise
    
    # Try to read version from file.
    if not found_git:
        try:
            if not os.path.exists('./data/'):
                os.makedirs('./data/')
            touch('./data/version')
            with open('./data/version', 'r') as f:
                version_by_git = f.readline().replace('\n', '')
                if version_by_git == '':
                    version_by_git = 'unkown'
        except:
            raise
    
    # Try to find sage version.
    if sage_path is not None:
        try:
            sage_version = subprocess.check_output( [ sage_path, '-v'], stderr=open(os.devnull, 'wb') ).replace('\n', ' ')
            found_sage = True
        except OSError:
            sage_version = "?"
    
    ret += 'Program version:  ' + version_by_git + '\n'
    ret += 'Operating System: ' + platform.system() + ' ' + platform.release() + '\n'
    ret += 'Python version:   ' + sys.version.replace('\n', ' ') + '\n'
    if sage_path is not None:
        ret += 'Sage version:     ' + sage_version + '\n'
    ret += 'Date:             ' + time.strftime("%Y-%b-%d %H:%M:%S") + '\n'
    
    return ret, True if found_sage or sage_path is None else False


def touch(filename, timestamp=None):
    """"Portable version of the GNU 'touch' command."""
    try:
        with open(filename, 'a'):
            os.utime(filename, timestamp)
    except:
        raise


def which(program):
    """"Portabel version of the GNU 'which' command."""
    def is_executable(path):
        return os.path.isfile(path) and os.access(path, os.X_OK)

    # Get Path and name of the program.
    # In case the path is empty, we search the progam using the PATH variable.
    path, name = os.path.split(program)
    if path:
        if is_executable(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            executable_file = os.path.join(path, program)
            if is_executable(executable_file):
                return executable_file
    return None


class Tee(object):
    """"A class simulating the linux program 'tee'."""
    def __init__(self, name, mode):
        self._file = None
        # Open file.
        try:
            self._file = open(name, mode)
        except:
            raise
        
        # Keep sys.stdout save.
        self._stdout = sys.stdout
        # Redirect sys.stdout to the instance of Tee.
        sys.stdout = self
    
    def __del__(self):
        self.stop_logging()
    
    def write(self, data):
        # Write to file.
        try:
            if self._file is not None:
                self._file.write(data)
        except:
            raise
        # Write to stdout.
        self._stdout.write(data)
    
    def flush(self):
        # Flush file.
        if self._file is not None:
            self._file.flush()
        # Flush stdout.
        self._stdout.flush()

    def stop_logging(self):
        # Restore sys.stdout.
        # It might happen, that sys is already shutdown.
        # Therefore, we try-catch this.
        try:
            sys.stdout = self._stdout
        except:
            pass
        # Close file.
        try:
            if self._file is not None:
                self._file.close()
                self._file = None
        except:
            pass


def write_chaincomplex_to_chomp_representation(open_file, diffs, verbose):
    """" This function is an alternation of the SageMath function
         def _chomp_repr_(self)" in "SageMath/local/lib/python2.7/site-packages/sage/homology/chain_complex.py"
         This provides an enormous speed up."""

    if len(diffs) == 0:
        diffs = {0: (0, 0, [])}
    maxdim = max(diffs)
    mindim = 0
    s = "chain complex\n\nmax dimension = %s\n\n" % (maxdim - mindim)
    try:
        open_file.write(s)
    except:
        raise

    for i in range(0, maxdim + 1):
        if verbose:
            sys.stdout.write("  Dimension " + str(i) + " ... ")
            sys.stdout.flush()
            starting_time = time.clock()

        s = "dimension %s\n" % i
        num_rows, num_cols, mat = diffs.get(i, (0, 0, []))

        # We sort the matrix dictionary which is  of the form
        #     (row, col) -> coeff
        # by the columns.
        sorted_keys = sorted(mat, key=lambda s: s[1])
        # col_idx is the column index treated last
        col_idx = -1
        # non_zero_column is set to True or Flase if the the treated column is zero or not
        non_zero_column = None

        # We iterate through columns via sorted keys
        for key in sorted_keys:
            # Check if the column index of the key changed
            if col_idx != key[1]:
                col_idx += 1
                # We finalize last column.
                # Observe that we use True and False
                # In particular we do not print an extra new line on the first iteration on sorted_keys.
                if non_zero_column is True:
                    s += "\n"
                elif non_zero_column is False:
                    s += "0\n"
                # We write empty columns (if any)
                while col_idx < key[1]:
                    s += "   boundary a%s = 0\n" % (col_idx + 1)
                    col_idx += 1
                # We start next column
                s += "   boundary a%s = " % (col_idx + 1)
                non_zero_column = False
            # Get the coefficient and write the entry if it is non-zero
            coeff = mat[key]
            if coeff > 0:
                s += "+ %s * a%s " % (coeff, key[0] + 1)
                non_zero_column = True
            elif coeff < 0:
                s += "- %s * a%s " % (-coeff, key[0] + 1)
                non_zero_column = True

        # We finalize last column that was treated by iteration through sorted_keys.
        # If sorted_keys was not empty then non_zero_column is either True or False. Otherwise it is None.
        if non_zero_column is True:
            s += "\n"
        elif non_zero_column is False:
            s += "0\n"
        col_idx += 1

        # We treat the remaining zero columns.
        while col_idx < num_cols:
            s += "   boundary a%s = 0\n" % (col_idx + 1)
            col_idx += 1
        s += "\n"

        try:
            open_file.write(s)
        except:
            raise

        if verbose:
            sys.stdout.write('Done. Duration = ' + str(time.clock() - starting_time) + '\n')
            sys.stdout.flush()


def write_performance_test(diffs):
    with open('performance_test.py', 'w') as datei:
        datei.write('from sage.all import ChainComplex, matrix, ZZ\n')
        datei.write('dict_chaincomplex = dict()\n')
        for degree, mat in diffs.items():
            num_rows, num_cols, bdry_matrix_dict = mat
            datei.write(
                'dict_chaincomplex[{degree}] = matrix(ZZ, {num_rows}, {num_cols}, {bdry_matrix_dict}, sparse=True)\n'.format(
                    degree=degree, num_rows=num_rows, num_cols=num_cols, bdry_matrix_dict=bdry_matrix_dict))
        datei.write('cplx = ChainComplex(dict_chaincomplex, degree_of_differential=-1, check=True)\n')
        datei.write('with open("test.bin", "w") as f:\n')
        datei.write('  f.write(cplx._chomp_repr_())\n')