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

from itertools import permutations
from gi.repository import Gtk, Gdk
from pyradbar.cell import *
from pyradbar.cell_para import *
from pyradbar.draw_cell import *
from pyradbar.perm import *
from cell_viewer import *

import os
import sys
import multiprocessing
import numpy.random
import gzip
import cPickle as pickle
import inspect
import itertools
import functools
import subprocess

def process_stuff(par):
    g = par[0]
    m = par[1]
    A = permutations( range(2*m-2+4*g) )
    create_directories(g, m)
    
    for a in A:
        cell = Cell(g, m, a)
        if( cell ):
            print( cell )
            #show_cell = ShowCell( cell )
            #Gtk.main()
            draw_cell = DrawCell( cell )

def parallelized_test():
    num_threads = 4
    pool = multiprocessing.Pool(num_threads)
    
    para_lst  = zip( 7*[0], [ x+1 for x in range(7) ] )
    para_lst += zip( 5*[1], [ x+1 for x in range(5) ] )
    para_lst += zip( 3*[2], [ x+1 for x in range(3) ] )
    para_lst += zip( 1*[3], [ x+1 for x in range(1) ] )
    
    pool.map(process_stuff, para_lst)

def face_map_test(cell):
    if cell is not None :
        for i in range( cell.degree() + 1 ):            
            bdry = cell.get_clean_copy()
            bdry.face(i)
            CellViewer( cell=bdry )
    else:
        while True :
            k = numpy.random.randint(low=2, high=12)
            tau = numpy.random.permutation(k)
            cell = Cell(tau)
            if cell.get_genus() == 1:
                print "Starting with new cell."
                while cell :
                    print( cell )
                    show_cell = DrawCell( cell )
                    Gtk.main()
                    i = numpy.random.randint(low=0, high=1+(k if k > 0 else 1))
                    k = k - 1
                    cell.face(i)
                print "Cell degenerate."

def store_taus(h):
        if h < 1:
            print "2*g+m is to small"
            return
        
        # m and archive as a function of g
        m_of = []
        archive_of = []
        for g in range( (h+2)/ 2 ):
            #exclude m = 0
            if h - 2*g <= 0:
                continue
            
            m_of.append( h-2*g )
            try:
                create_directories(g, m_of[g])
                archive_of.append( gzip.GzipFile('./data/top_cell_g_' + str(g) + '_m_' + str(m_of[g]) + '.bz2', 'wb') )
            except:
                # Print the error.
                frameinfo = inspect.getframeinfo(inspect.currentframe())
                print frameinfo.filename, frameinfo.lineno
                e, p, t = sys.exc_info()
                print e, p
                return
        
        # iterate through permutations
        taus = permutations( range(2*h-2) )
        total = math.factorial(2*h-2)
        cnt = 0
        for tau in taus:
            cnt += 1
            if cnt % 10000 == 0:
                sys.stdout.write( "\r" + str(100*cnt/total) + "%" )
                sys.stdout.flush()
            cell = Cell(tau)
            if cell.is_top_cell() == True:
                g = cell.get_genus()
                m = m_of[g]
                pickle.dump( tau, archive_of[g] )
        sys.stdout.write("\r")

def hash_tau( tau, dict=None ):
    if dict is None:
        dict = {}
    cell = Cell(tau)
    dict[cell] = dict.get(cell, len(dict))

def cache_cells():
    for h in range(8):
        print "Processing h =", h
        store_taus(h)

def test_loading(g=1, m=2):
    local_dict = {}
    with LoadTopCellRho(g, m) as rho_archive:
        for rho in rho_archive:
            hash_tau(rho, local_dict)
    
    print "g =", g, "m =", m, "num cells =", len(local_dict)
    
    for key,val in local_dict.items():
        print key, "->", val

def test_loading_all():
    for h in range(8):
        for g in range( (h+2)/ 2 ):
            #exclude m = 0
            if h - 2*g <= 0:
                continue
            
            test_loading(g, h-2*g)

def show_faces(g=0, m=2):
    print "Loading top dimensional cells."
    next_basis = {}
    dict_chaincomplex = {}
    degree = 4*g+2*m
    with LoadTopCellTau(g, m) as tau_archive:
        for tau in tau_archive:
            cell = Cell(tau)
            print cell
            next_basis[cell] =  next_basis.get(cell, len(next_basis))
    print " "
    
    next_basis, num_rows, num_cols = compute_faces_and_draw(next_basis)
    for cell in next_basis:
        print cell
    print " "
    while (next_basis is not None) and (len(next_basis) > 0):
        next_basis, num_rows, num_cols = compute_faces_and_draw(next_basis)
        for cell in next_basis:
            print cell
        print " "
        degree -= 1

def compute_faces_and_draw( cells, print_coefficients=False ):
    next_basis = {}
    num_rows = 0
    num_cols = 0
    
    if cells is None or len(cells) == 0:
        return next_basis, num_rows, num_cols
    
    # The 'first' element of a dictionary is given by cells.iterkeys().next()
    degree = cells.iterkeys().next().degree()
    if degree == 0:
        return next_basis, num_rows, num_cols
    
    for cell in cells:
        cell_idx = cells[cell]
        # compute the boundary of the given cell
        coefficients = {}
        sign = 1;
        #ShowCell(cell)
        #Gtk.main()
        for i in range( degree+1 ):
            face = cell.get_clean_copy()
            face.face(i)
            #print "    " + str(face)
            #ShowCell(face)
            #Gtk.main()
            face_idx = next_basis[face] = next_basis.get(face, len(next_basis))
            coefficients[face_idx] = coefficients.get( face_idx, 0 ) + sign
            sign *= -1
        if print_coefficients:
            for key, val in enumerate(coefficients):
                print key, '->', val
    
    num_cols = len(cells)
    num_rows = len(next_basis)
    
    return next_basis, num_rows, num_cols

def test_cell_para():
    rho = (3, 2, 1, 0)
    marking = [ [0, 1], [2, 2] ]
    cell = CellPara( rho, marking )
    print cell
    cell.face(3)
    print cell
    cell.face(1)
    print cell
    cell.face(1)
    print cell

def test_cell_parametrizations():
    
    rho = (2, 1, 0)
    cell = Cell(rho, suspend=False)
    m = cell.get_num_punc()
    print cell
    
    for position_tuple in list( itertools.product( *cell.cyc_rho() ) ):
        positions = list(position_tuple)
        positions.sort()
        marks = m*[0]
        offset = 0
        p = len(cell.get_rho()) + len(positions)
        rho = p*[0]
        for i, x in enumerate(cell.get_rho()):
            rho[ i + offset ] = x + len( [ a for a in positions if x >= a ] )
            if i in positions:
                offset += 1
                rho[ i + offset ] = i + offset - 1
                marks[offset - 1] = i + offset
        print positions, marks, rho
        print ""
        
        #for marking_perm in permutations( range(m) ):
        #    marking = { marks[x] : 1 + marking_perm[x] for x in range(m) }
        #    print CellPara( rho, marking )

def test_cell_parametrizations_factory():
    rho = (0, 3, 2, 1)
    unparametrized_cell = Cell(rho, suspend=False)
    for cell in ParametrizeBoundary(unparametrized_cell):
        print cell

def test_noncrossing_partitions():
    dw = DyckWords(6)
    for w in dw:
        p = w.to_noncrossing_partition()
        valid=True
        for l in p:
            if len(l) != 2:
                valid=False
                break
        if valid:
            print p

def main():
    #test_cell_parametrizations_factory()
    #test_noncrossing_partitions()
    #test_loading(g=0, m=6)
    face_map_test( Cell( (0, 6, 4, 5, 3, 1, 2), suspend=False ) )
    
if __name__ == "__main__":
    main()
