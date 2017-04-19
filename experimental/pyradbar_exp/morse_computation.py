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
import collections
import copy
import sys

def noncrossing_partitions_of_tuples(number_of_symbols=None):
    lst_noncrossing_partitions_of_tuples = []
    # There are not partitions of this type for an odd number of symbols.
    if number_of_symbols % 2 == 1:
        return lst_noncrossing_partitions_of_tuples
    
    # We generate all noncrossing partitions from the corresponding Dyck words.
    dw = DyckWords(number_of_symbols)
    for w in dw:
        noncrossing_partition = w.to_noncrossing_partition()
        
        # We are interested only in those consisting of blocks of size 2.
        consists_of_tuples=True
        for block in noncrossing_partition:
            if len(block) != 2:
                consists_of_tuples=False
                break
        if consists_of_tuples == True:
            # Renormalize.
            valid_partition = []
            for block in noncrossing_partition:
                valid_partition += [[ block[0] - 1,  block[1] - 1 ]]
            lst_noncrossing_partitions_of_tuples += [valid_partition]
    
    return lst_noncrossing_partitions_of_tuples

def compute_faces_matrix_with_morse_order( old_basis, degree ):
    #
    # Compare the compute_face_matrix function from homology_computation.py
    #
    
    # Setup variables.
    next_basis = { 'e' : {}, 'c' : {}, 'r' : {} }
    num_rows = 0
    num_cols = 0
    
    # Return empty stuff if input is empty.
    if old_basis is None or degree == 0:
        return next_basis, num_rows, num_cols
    
    # Compute boundaries cell by cell.
    cells = [] 
    for t in {'e', 'c', 'r'}:
        cells.extend( old_basis[t].keys() )
    
    for cell in cells:
        for i in range( degree+1 ):
            face = cell.get_clean_copy()
            face.face(i)
            mt, j = face.morse_type()
            next_basis[mt][face] = next_basis[mt].get(face, len(next_basis[mt]))
    
    # Done.
    return next_basis, num_rows, num_cols

def print_morse_differential(essentials_high, essentials_low, graph):
    num_cols = len(essentials_high)
    num_rows = len(essentials_low)
    matrix_dict = {}
    
    if num_cols == 0 or num_rows == 0:
        return
    
    graph.add_vertices(essentials_high)
    graph.add_vertices(essentials_low)
    
    for cell_high, j in essentials_high.items():
        for cell_low, i in essentials_low.items():
            coefficient = 0
            for path in graph.all_paths( cell_high, cell_low ):
                path_coefficient = 1
                for k in range(len(path)-1):
                    path_coefficient *= graph.edge_label(path[k], path[k+1])
                coefficient += path_coefficient
            matrix_dict[i, j] = coefficient
    
    morse_differential = matrix( ZZ, num_rows, num_cols, matrix_dict, sparse=True )
    
    sys.stdout.write( morse_differential.str().replace(' 0', ' .').replace('[0', '[.') + '\n')
    

def show_morse_order( old_basis, next_basis, degree, print_essentials=False ):
    essentials = old_basis['e']
    collapsibles = old_basis['c']
    redundants = next_basis['r']
    
    # Print the degree.
    sys.stdout.write('\nDegree = ' + str(degree) + '\n')
    
    # Print the essential cells if wanted.
    if len(essentials) > 0:
        sys.stdout.write('Number of essentials = ' + str(len(essentials)) + '\n')
        if print_essentials:
            for cell in essentials.keys():
                sys.stdout.write(str(cell)+'\n')
    
    # Create the graph.
    graph = DiGraph(weighted=True)
    
    # Prepare the weighted edges.
    edge_list = []
    
    # We add the collapsibles and also their redundant and essential faces.
    for cell in collapsibles.keys():
        t, j = cell.morse_type()
        
        # Compute boundary.
        coefficients = {}
        sign = 1;
        for i in range( degree+1 ):
            face = cell.get_clean_copy()
            face.face(i)
            # it suffices to consider redundant cells.
            mt, midx = face.morse_type()
            if mt in {'e', 'r'} :
                coefficients[face] = coefficients.get( face, 0 ) + sign
            sign *= -1
        
        redundant_partner = cell.get_clean_copy()
        redundant_partner.face(j)
        
        # Add edges to the graph.
        for face, coeff in coefficients.items():
            if coeff != 0:
                if face == redundant_partner:
                    edge_list.append( (redundant_partner, cell, -coeff ) )
                else:
                    edge_list.append( (cell, face, coeff ) )
    
    # We add the essentials and also their redundant and essential faces.
    for cell in essentials.keys():
        # Compute boundary.
        coefficients = {}
        sign = 1;
        for i in range( degree+1 ):
            face = cell.get_clean_copy()
            face.face(i)
            # it suffices to consider redundant cells.
            mt, midx = face.morse_type()
            if mt in {'e', 'r'} :
                coefficients[face] = coefficients.get( face, 0 ) + sign
            sign *= -1
        
        # Add edges to the graph.
        for face, coeff in coefficients.items():
            if coeff != 0:
                edge_list.append( (cell, face, coeff ) )
    
    # Add the weighted edges to the graph.
    graph.add_edges(edge_list)
    sys.stdout.write(str(graph) + '\n')
    
    # Check the acyclicity condition.
    cycles_in_the_morse_graph = graph.all_simple_cycles()
    if cycles_in_the_morse_graph is not None and len(cycles_in_the_morse_graph) > 0:
        sys.stdout.write( 'The Morse Graph has oriented cycles using cells in degrees ' + str(degree) + ' and ' + str (degree-1) + '.\n' )
        for cycle in cycles_in_the_morse_graph:
            for cell in cycle:
                sys.stdout.write( str(cell) + str(cell.morse_type()) + '\n' )
            sys.stdout.write( '\n-------------------+-------------------\n' )
        sys.stdout.write( '\n#############################################\n' )
        return
    else:
        sys.stdout.write( 'The Morse Graph has no oriented cycles using cells in degrees ' + str(degree) + ' and ' + str (degree-1) + '.\n' )
    
    # Print the Morse differential.
    print_morse_differential(old_basis['e'], next_basis['e'], graph)
    
    sys.stdout.write( '\n#############################################\n' )

def morse_test(g=0, m=5, more_verbose=False):
    more_verbose = True if more_verbose == 'True' else False
    next_basis = { 'e' : {}, 'c' : {}, 'r' : {} }
    degree = 4*g+2*m-1
    
    if g == 0 and m > 7:
        p = 4*g+2*m-2
        partitions = noncrossing_partitions_of_tuples(2*m-2)
        for part in partitions:
            lam = {}
            for part_tuple in part:
                lam[ part_tuple[0] ] = part_tuple[1]
                lam[ part_tuple[1] ] = part_tuple[0]
            rho = tuple( (lam[x] - 1 + p) % p for x in range(p) )
            
            cell = Cell(rho)
            mt, i = cell.morse_type()
            next_basis[mt][cell] = next_basis[mt].get(cell, len(next_basis[mt]))
    else:
        with LoadTopCellRho(g, m) as rho_archive:
            for rho in rho_archive:
                cell = Cell(rho)
                mt, i = cell.morse_type()
                next_basis[mt][cell] = next_basis[mt].get(cell, len(next_basis[mt]))
    degree -= 1
    
    do_while = True
    while do_while:
        old_basis = copy.deepcopy(next_basis)
        next_basis, num_rows, num_cols = compute_faces_matrix_with_morse_order(next_basis, degree)
        show_morse_order(old_basis, next_basis, degree, print_essentials=more_verbose)
        
        degree -= 1
        if degree < 0:
            do_while = False

def main(g=0, m=5, more_verbose=False, result_file=None):
    # Tee the output.
    tee = Tee(result_file, 'a' )
    
    if g > 0:
        sys.stdout.write('\nFor g > 0, the algorithm is not yet working correctly.\n\n')
    
    morse_test(g, m, more_verbose)
    sys.stdout.flush()

main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4])
