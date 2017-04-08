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

# Print a nice preamble.
def preamble( sage_path = None ):
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
            with open('./data/version', 'rb') as f:
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

# Portable version of the GNU 'touch' command.
def touch(filename, timestamp=None):
    try:
        with open(filename, 'a'):
            os.utime(filename, timestamp)
    except:
        raise

# Portabel version of the GNU 'which' command.
def which(program):
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

# A class simulating the linux program 'tee'.
class Tee(object):
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
