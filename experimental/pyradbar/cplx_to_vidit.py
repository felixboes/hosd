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

# Note: The script has to be called from sage itself
# sage -python /path/to/script.py
from sage.all import *
from cell import *
from misc_useability_stuff import *
import sys
import time
import collections 
import copy

def cplx_to_vidit(vidit_file_prefix, g = 1, m = 2, ring='ZZ', verbose=True):
    # Setup coefficients and other variables.
    verbose = True if verbose == 'True' else False
    coeff_ring = None
    try:
        coeff_ring = eval(ring)
    except:
        sys.stdout.write('\n')
        sys.stdout.write('    Error: Could not identify your ring \''+str(ring)+'\'. Defaulting to the rationals QQ.\n')
        sys.stdout.write('\n')
        coeff_ring = QQ
    
    # Setup other variables.
    next_basis = {}
    dict_chaincomplex = {}
    degree = 4*g + 2*m-2
    starting_time = None
    
    # Load top cells.
    if verbose:
        sys.stdout.write('We construct the cellular complex for g = ' + str(g) + ' and m = ' + str(m) + ' with coefficients in \''+str(coeff_ring)+'\'.\n')
        sys.stdout.write('Then we save the complex to file.\n\n')
        sys.stdout.write('Loading top cells ... ')
        sys.stdout.flush()
        starting_time = time.clock()
    with LoadTopCellRho(g, m) as rho_archive:
        for rho in rho_archive:
            cell = Cell(rho)
            next_basis[cell] =  next_basis.get(cell, len(next_basis))
    if verbose:
        sys.stdout.write('Done. Duration = ' + str(time.clock() - starting_time) + '\n')
    
    # Compute the differentials and store them to file (one file per matrix).
    while degree >= 0:
        # Open file.
        filename = vidit_file_prefix + '_' + str(degree)
        if verbose:
            sys.stdout.write('Opening file \'' + filename + '\' ... ')
            sys.stdout.flush()
        try:
            open_file = open(filename, 'wb')
        except:
            sys.stdout.write(' ABROATING! There was an error while opening the file ' + filename + '.\n')
            frameinfo = inspect.getframeinfo(inspect.currentframe())
            e, p, t = sys.exc_info()
            sys.stdout.write( str(frameinfo.filename) + ' ' + str(frameinfo.lineno) + '\n' + str(e) + ' ' + str(p) + '\n')
            sys.stdout.flush()
            raise
        else:
            sys.stdout.write('Done.\n')
            sys.stdout.flush()
        
        # Compute Differential and save to file.
        if verbose:
            sys.stdout.write('Computing the differential D_' + str(degree) + ' ... ')
            sys.stdout.flush()
            starting_time = time.clock()
        try:
            next_basis = compute_faces_for_vidit(next_basis, open_file)
            open_file.close()
        except:
            sys.stdout.write(' ABROATING! There was an error while writing the file ' + filename + '.\n')
            frameinfo = inspect.getframeinfo(inspect.currentframe())
            e, p, t = sys.exc_info()
            sys.stdout.write( str(frameinfo.filename) + ' ' + str(frameinfo.lineno) + '\n' + str(e) + ' ' + str(p) + '\n')
            sys.stdout.flush()
            try:
                open_file.close()
            except:
                pass
            raise
        degree -= 1
        if verbose:
            sys.stdout.write('Done. Duration = ' + str(time.clock() - starting_time) + '\n')
    
    # Concaternate all files.
    degree = 0
    if verbose:
        sys.stdout.write('Concaternating files.')
        sys.stdout.flush()
    try:
        filename = vidit_file_prefix
        out_file = open(filename, 'wb')
    except:
        sys.stdout.write(' ABROATING! There was an error while opening the file ' + filename + '.\n')
        frameinfo = inspect.getframeinfo(inspect.currentframe())
        e, p, t = sys.exc_info()
        sys.stdout.write( str(frameinfo.filename) + ' ' + str(frameinfo.lineno) + '\n' + str(e) + ' ' + str(p) + '\n')
        sys.stdout.flush()
        raise
    else:
        sys.stdout.write('Done.\n')
        sys.stdout.flush()
    
    input_filenames = [ vidit_file_prefix + '_' + str(d) for d in range(4*g + 2*m - 1) ]
    try:
        for filename in input_filenames:
            with open(filename) as in_file:
                for line in in_file:
                    out_file.write(line)
            try:
                os.remove(filename)
            except:
                sys.stdout.write(' ABROATING! There was an error while removing ' + filename + '.\n')
                frameinfo = inspect.getframeinfo(inspect.currentframe())
                e, p, t = sys.exc_info()
                sys.stdout.write( str(frameinfo.filename) + ' ' + str(frameinfo.lineno) + '\n' + str(e) + ' ' + str(p) + '\n')
                sys.stdout.flush()
                raise
        out_file.close()
    except:
        sys.stdout.write(' ABROATING! There was an error while opening / reading / writing files around ' + filename + '.\n')
        frameinfo = inspect.getframeinfo(inspect.currentframe())
        e, p, t = sys.exc_info()
        sys.stdout.write( str(frameinfo.filename) + ' ' + str(frameinfo.lineno) + '\n' + str(e) + ' ' + str(p) + '\n')
        sys.stdout.flush()
        raise

def compute_faces_for_vidit( cells, open_file ):
    # For Vidits program, we have to provide a homological chain complex from botton to top.
    # If want give a filtration, we have to indicate the filtration steps.
    # Right now, we don't do it.
    #
    # The line 
    # -1
    # indicates that we work with the next dimension.
    #
    # In each dimension, every line describes a unique basis element and the linear combination of its faces.
    # 1 k a_1 f_1 a_2 f_2 ... a_k f_k
    # where
    #   1 = filtration index (we do not care about right now).
    #   k = the length of the linear combination.
    #   a_i = the coefficient.
    #   f_i = the f_i th basis element which is hit.
    # The enumeration of each basis starts with zero.
    #
    # Example: We provide the (standard) cell complex for RP^2, with 1 zero, one and two cell.
    # -1
    # 1 0
    # -1
    # 1 0
    # -1
    # 1 1 2 1
    
    # Setup variables.
    next_basis = {}
    
    # New basis.
    open_file.write('-1\n')
    
    # Return empty stuff if input is empty.
    if cells is None or len(cells) == 0:
        return next_basis
    
    # Get the degree.
    # The 'first' element of a dictionary is given by cells.iterkeys().next()
    degree = cells.iterkeys().next().degree()
    if degree == 0:
        # There is a single zero cell.
        open_file.write('1 0\n')
        return next_basis
    
    # Compute boundaries cell by cell.
    for cell in cells:
        cell_idx = cells[cell]
        # Compute the boundary of the given cell.
        coefficients = {}
        sign = 1;
        for i in range( degree+1 ):
            face = cell.get_clean_copy()
            face.face(i)
            face_idx = next_basis[face] = next_basis.get(face, len(next_basis))
            coefficients[face_idx] = coefficients.get( face_idx, 0 ) + sign
            sign *= -1
        
        # Prepare linear combination.
        num_bdry_cells = 0
        lin_comb = ''
        for face_idx, coeff in coefficients.items():
            if coeff != 0:
                num_bdry_cells += 1
                lin_comb += ' ' + str(coeff) + ' ' + str(face_idx)
        
        # Store linear combination.
        open_file.write('1 ' + str(num_bdry_cells) + lin_comb + '\n')
    
    # Done.
    return next_basis

def main(vidit_file_prefix, g=0, m=7, ring=None, verbose=None ):
    cplx_to_vidit(vidit_file_prefix, g, m, ring, verbose )

main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5] )
