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

from cell import *
from cell_para import *
from misc_useability_stuff import *
import sys
import time

def compute_homology(g = 1, m = 2, param_bdry=False, verbose=True, more_verbose=False, only_good_stuff=False, homchain_file=None ):
    # Setup coefficients and other variables.
    if g < 0 or m <= 0:
        sys.stdout.write('The genus has to be non-negative and the number of punctures has to be positive.\nGot g = ' + str(g) + '; m = ' + str(m) + '.\n')
        sys.stdout.write('ABROATING!\n')
        sys.stdout.flush()
        
        raise ValueError()
    
    param_bdry = True if param_bdry == 'True' else False
    verbose = True if verbose == 'True' else False
    more_verbose = True if more_verbose == 'True' else False
    only_good_stuff = True if only_good_stuff == 'True' else False
    open_file = None
    sys.stdout.write('Opening file \'' + homchain_file + '\' ... ')
    sys.stdout.flush()
    try:
        open_file = open(homchain_file, 'wb')
    except:
        sys.stdout.write(' ABROATING! There was an error while writing the homchain file ' + homchain_file + '.\n')
        frameinfo = inspect.getframeinfo(inspect.currentframe())
        e, p, t = sys.exc_info()
        sys.stdout.write( str(frameinfo.filename) + ' ' + str(frameinfo.lineno) + '\n' + str(e) + ' ' + str(p) + '\n')
        sys.stdout.flush()
    else:
        sys.stdout.write('Done.\n')
        sys.stdout.flush()
    
    # Setup other variables.
    next_basis = {}
    dict_chaincomplex = {}
    degree = 4*g + 2*(m-1) if param_bdry == False else 4*g + 3*(m-1) + 1
    starting_time = None
    
    # Load top cells.
    if verbose:
        sys.stdout.write('We construct the cellular complex for g = ' + str(g) + ' and m = ' + str(m) + ' with ' + ( 'un' if param_bdry == False else '')  +'parametrized boundaries.\n')
        sys.stdout.write('Then we save the chain complex in chomp representation to \'' + homchain_file + '\'.\n\n')
        sys.stdout.write('Loading top cells ... ')
        sys.stdout.flush()
        starting_time = time.clock()
    if g > 0 or m > 0:
        LoadTopCells = (LoadTopCellRho(g,m) if param_bdry == False else LoadTopCellParaRho(g,m))
        with LoadTopCells as rho_archive:
            # Check wether the archive exists or not.
            if rho_archive is None:
                sys.stdout.write(
                    "Loading top cells failed. We assume you did not create run our program 'create_and_store_top_cells.py' that creates and stores top cells for a given number hh = 2g+m.\n")
                sys.stdout.write('\n\n\n      ABROATING!\n\n\n')
                sys.stdout.flush()
                raise RuntimeError()

            for rho in rho_archive:
                cell = (Cell(rho, suspend=True) if param_bdry == False else CellPara(rho, degree - 1, suspend=True))
                next_basis[cell] = next_basis.get(cell, len(next_basis))
    if verbose:
        sys.stdout.write('Done. Duration = ' + str(time.clock() - starting_time) + '\n')
    
    # Compute the first differential.
    if verbose:
        sys.stdout.write('Computing the differential D_' + str(degree) + ' ... ')
        sys.stdout.flush()
        starting_time = time.clock()

    next_basis, num_rows, num_cols, bdry_matrix_dict = compute_faces_matrix(next_basis, only_good_stuff)
    dict_chaincomplex[degree] = (num_rows, num_cols, bdry_matrix_dict)

    degree -= 1

    if verbose:
        sys.stdout.write('Done. Duration = ' + str(time.clock() - starting_time) + '\n')
    
    # Compute the remaining differentials.
    while (next_basis is not None) and (len(next_basis) > 0):
        if verbose:
            sys.stdout.write('Computing the differential D_' + str(degree) + ' ... ')
            sys.stdout.flush()
            starting_time = time.clock()

            next_basis, num_rows, num_cols, bdry_matrix_dict = compute_faces_matrix(next_basis, only_good_stuff)
            dict_chaincomplex[degree] = (num_rows, num_cols, bdry_matrix_dict)

        if verbose:
            sys.stdout.write('Done. Duration = ' + str(time.clock() - starting_time) + '\n')
        degree -= 1
    
    # Construct the chain complex from the dictionaries.
    if verbose:
        sys.stdout.write('Constructing the chain complex ... ')
        sys.stdout.flush()
        starting_time = time.clock()

    sys.stdout.write('Writing chain complex representation to file \'' + homchain_file + '\' : \n')
    sys.stdout.flush()
    try:
        write_chaincomplex_to_chomp_representation( open_file, dict_chaincomplex, verbose )
    except:
        sys.stdout.write(' ABROATING! There was an error while writing the homchain file ' + homchain_file + '.\n')
        frameinfo = inspect.getframeinfo(inspect.currentframe())
        e, p, t = sys.exc_info()
        sys.stdout.write( str(frameinfo.filename) + ' ' + str(frameinfo.lineno) + '\n' + str(e) + ' ' + str(p) + '\n')
        sys.stdout.flush()
        try:
            open_file.close()
        except:
            pass
        raise
    else:
        try:
            open_file.close()
        except:
            raise
    return

def compute_faces_matrix( cells, only_good_stuff=False ):
    # Recall the following passage from the sage 6.7 manual:
    # The entries of a matrix can be specified as a flat list of elements,
    # a list of lists (i.e., a list of rows), a list of Sage vectors, a callable object,
    # or a dictionary having positions as keys and matrix entries as values (see the examples).
    #
    # We are going to use a dictionary.
    
    # Setup variables.
    next_basis = {}
    matrix_dict = {}
    num_rows = 0
    num_cols = 0
    
    # Return empty stuff if input is empty.
    if cells is None or len(cells) == 0:
        return next_basis, num_rows, num_cols, matrix_dict
    
    # Get the degree.
    # The 'first' element of a dictionary is given by cells.iterkeys().next()
    degree = cells.iterkeys().next().degree()
    if degree == 0:
        num_cols = len(cells)
        return next_basis, num_rows, num_cols, matrix_dict
    
    # Compute boundaries cell by cell.
    for cell in cells:
        cell_idx = cells[cell]
        # compute the boundary of the given cell
        coefficients = {}
        sign = 1
        for i in range( degree+1 ):
            face = cell.get_clean_copy()
            face.face(i)
            if only_good_stuff == True:
                if face.get_degenerated() == False:
                    face_idx = next_basis[face] = next_basis.get(face, len(next_basis))
                    coefficients[face_idx] = coefficients.get( face_idx, 0 ) + sign*cell.orientation_sign(i)
            else:
                face_idx = next_basis[face] = next_basis.get(face, len(next_basis))
                coefficients[face_idx] = coefficients.get( face_idx, 0 ) + sign
            sign *= -1
        
        # store the column in the dictionary
        for face_idx, coeff in coefficients.items():
            matrix_dict[ face_idx, cell_idx ] = coeff
    
    # Compute number of columns and rows.
    num_cols = len(cells)
    num_rows = len(next_basis)

    # Done.
    return next_basis, num_rows, num_cols, matrix_dict

def main(g=0, m=7, param_bdry=False, verbose=None, more_verbose=None, only_good_stuff=None, result_file=None, homchain_file=None):
    # Tee the output.
    tee = Tee(result_file, 'a' )
    
    compute_homology(g, m, param_bdry, verbose, more_verbose, only_good_stuff, homchain_file)
    sys.stdout.write('The Homology of the complex should be computed using the program \'homchain_gmp\'.\n')
    sys.stdout.flush()

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
