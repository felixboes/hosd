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

import inspect
import sys
import time

from .symmetric_group import Permutation, SymmetricGroup
from .misc_useability_stuff import write_chaincomplex_to_chomp_representation

def compute_homology_symm(max_d, verbose=True, homchain_file=None):
  # Open homchain file.
  open_file = None
  if homchain_file is not None:
    print('Opening file \'' + homchain_file + '\' ... ', end='')
    try:
      open_file = open(homchain_file, 'w')
    except:
      print(' ABROATING! There was an error while writing the homchain file', homchain_file)
      frameinfo = inspect.getframeinfo(inspect.currentframe())
      e, p, t = sys.exc_info()
      print(str(frameinfo.filename), str(frameinfo.lineno))
      print(str(e), str(p))
    else:
      print('Done.')

    # Setup other variables.
    next_basis = {}
    dict_chaincomplex = {}
    degree = max_d
    starting_time = None

    # Load top cells.
    if verbose:
      print('We construct the cellular complex for Symm_d with d <= ', degree, '.', sep='')
      print('Then we save the chain complex in chomp representation to ', homchain_file, '.', sep='')
      print('')
      print('Loading top cells ... ', end='')
      starting_time = time.process_time()

    sym = SymmetricGroup(max_d)
    cells = sym.get_long_cycles()
    for cell in cells:
      next_basis[cell] = next_basis.get(cell, len(next_basis))

    if verbose:
      print('Done. Duration = ', time.clock() - starting_time, sep='')

    # Compute the first differential.
    if verbose:
      print('Computing the differential D_', degree, ' ... ', sep='', end='')
      starting_time = time.clock()

    next_basis, num_rows, num_cols, bdry_matrix_dict = compute_faces_matrix(next_basis, degree)
    dict_chaincomplex[degree] = (num_rows, num_cols, bdry_matrix_dict)

    degree -= 1

    if verbose:
      print('Done. Duration = ', time.clock() - starting_time, sep='')

    # Compute the remaining differentials.
    while (next_basis is not None) and (len(next_basis) > 0):
      if verbose:
        print('Computing the differential D_', degree,' ... ', sep='', end='')
        starting_time = time.clock()

      next_basis, num_rows, num_cols, bdry_matrix_dict = compute_faces_matrix(next_basis, degree)
      dict_chaincomplex[degree] = (num_rows, num_cols, bdry_matrix_dict)

      if verbose:
        print('Done. Duration = ', time.clock() - starting_time, sep='')
      degree -= 1

    # Construct the chain complex from the dictionaries.
    print('Writing chain complex representation to file ', homchain_file, ':', sep='')
    try:
      write_chaincomplex_to_chomp_representation(open_file, dict_chaincomplex, verbose)
    except:
      print(' ABROATING! There was an error while writing the homchain file', homchain_file)
      frameinfo = inspect.getframeinfo(inspect.currentframe())
      e, p, t = sys.exc_info()
      print(str(frameinfo.filename), str(frameinfo.lineno))
      print(str(e), str(p))
    finally:
      open_file.close()


def compute_faces_matrix(cells, degree):
  # Setup variables.
  next_basis = {}
  matrix_dict = {}
  num_rows = 0
  num_cols = 0

  # Return empty stuff if input is empty.
  if cells is None or len(cells) == 0:
    return next_basis, num_rows, num_cols, matrix_dict

  if degree == 0:
    num_cols = len(cells)
    return next_basis, num_rows, num_cols, matrix_dict

  # Compute boundaries cell by cell.
  for cell in cells:
    cell_idx = cells[cell]
    # compute the boundary of the given cell
    coefficients = {}
    sign = 1
    for i in range(degree + 1):
      face = cell.face(i)
      face_idx = next_basis[face] = next_basis.get(face, len(next_basis))
      coefficients[face_idx] = coefficients.get(face_idx, 0) + sign
      sign *= -1

    # store the column in the dictionary
    for face_idx, coeff in coefficients.items():
      matrix_dict[face_idx, cell_idx] = coeff

  # Compute number of columns and rows.
  num_cols = len(cells)
  num_rows = len(next_basis)

  # Done.
  return next_basis, num_rows, num_cols, matrix_dict
