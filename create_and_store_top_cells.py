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
import cPickle as pickle
import gzip
import inspect
import sys
import time

import pyradbar

# Produce a list of all indices in T_1 x T_2 x ... T_m
def count( list_of_tuples ):
    if list_of_tuples is not None and len(list_of_tuples) == 0:
        return []
    if len(list_of_tuples) == 1:
        return [ [i] for i in list_of_tuples[0] ]
    
    ret = []
    for i in list_of_tuples[0]:
        for remaining in count( list_of_tuples[1:] ):
            if remaining is not None:
                ret.append([i] + remaining)
    return ret

def parametrize( cyc_rhoprime, p, marking ):
    l_i = p + len(cyc_rhoprime)
    cyc_rho = [ [] for i in cyc_rhoprime ]
    for i, cyc_rhoprime_i in enumerate(cyc_rhoprime):
        for k in cyc_rhoprime_i:
            # Count the markings in front of k
            num_markings_in_front_of_k = len( [m for m in marking if k > m] )
            if k != marking[i]:
                cyc_rho[i].append( k + num_markings_in_front_of_k )
            else:
                cyc_rho[i].append( k + num_markings_in_front_of_k + 1 )
                cyc_rho[i].append( l_i )
                l_i += 1
                cyc_rho[i].append( k + num_markings_in_front_of_k )
    
    rho = (p + 2*len(cyc_rhoprime))*[0]
    for c in cyc_rho:
        for i in range(len(c)):
            rho[ c[i] ] = c[ (i + 1) % len(c) ]
    
    return rho

# Store all top cells for which 2g+m = h.
def store_top_cells( h, parametrized=False ):
    if h < 1:
        print "h = 2g+m-1 is to small."
        return
    
    sys.stdout.write('Creating and storing cells for 2g+m=' + str(h) + ' ... ' + '\r' )
    sys.stdout.flush()
    starting_time = time.clock()
    
    # m and archive as a function of g
    m_of = []
    archive_of = []
    for g in range( (h+2)/ 2 ):
        #exclude m = 0
        if h - 2*g <= 0:
            continue
        m_of.append( h-2*g )
        
        # Create directories.
        try:
            pyradbar.create_directories(g, m_of[g])
            archive_of.append( gzip.GzipFile('./data/top_cell_g_' + str(g) + '_m_' + str(m_of[g]) + ('_parametrized' if parametrized else '_unparametrized') + '.bz2', 'wb') )
        except:
            # Print the error.
            frameinfo = inspect.getframeinfo(inspect.currentframe())
            print frameinfo.filename, frameinfo.lineno
            e, p, t = sys.exc_info()
            print e, p
            return
    
    # iterate through all pairs
    def list_of_pairs( lst ):
        if len(lst) < 2:
            yield lst
            return
        low = lst[0]
        for i in range(1,len(lst)):
            pair = (low,lst[i])
            for remaining in list_of_pairs( lst[1:i] + lst[i+1: ] ):
                yield [pair] + remaining
    
    lams = list( list_of_pairs( range(0, 2*h-2) ) )
    l = len(lams)
    for n, lam in enumerate(lams):
        lamprime = (2*h-2)*[0]
        for t in lam:
            lamprime[ t[0] ] = t[1]
            lamprime[ t[1] ] = t[0]
        rhoprime = tuple( (lamprime[x] + 2*h - 3) % (2*h-2) for x in range(2*h-2) )
        
        if n % 5000 == 0:
            sys.stdout.write('Creating and storing cells for 2g+m=' + str(h) + ' ... '  + "{:.0%}".format(float(n)/l) + '\r')
            sys.stdout.flush()
        
        if parametrized == False:
            rho = rhoprime
            # Create Cell from permutation. This cell is a top cell by construction.
            cell = pyradbar.Cell(rho)
            g = cell.get_genus()
            pickle.dump( rho, archive_of[g] )
        else:
            cyc_rhoprime = pyradbar.cycle_decomposition(rhoprime)
            p = len(rhoprime)
            for i in count( cyc_rhoprime ):
                # Generate rho from rho' and a choice of markers in each outgoing boundary.
                rho = parametrize( cyc_rhoprime, p, i )
                deg = p + len(cyc_rhoprime) - 1
                # Create Cell from permutation. This cell is a top cell by construction.
                cell = pyradbar.CellPara( rho,  deg )
                g = cell.get_genus()
                pickle.dump(rho, archive_of[g])

                # Rotate counter clock wise, if there is a marking at p-1.
                if p-1 in i:
                    # shift suspension to 0.
                    r = len(rho)*[0]
                    for x in range(deg+1):
                        symbol = rho[ (x+(deg+1)-1)%(deg+1) ]
                        r[x] = symbol 
                        if symbol < (deg+1):
                            r[x] = (r[x]+1)% (deg+1)
                    for x in range(deg+1, len(rho)):
                        r[x] = (rho[x] + 1 ) % (deg+1)
                    
                    rho = r
                    # Create Cell from permutation. This cell is a top cell by construction.
                    cell = pyradbar.CellPara(rho, deg)
                    g = cell.get_genus()
                    pickle.dump(rho, archive_of[g])
    
    sys.stdout.write('Creating and storing cells for 2g+m=' + str(h) + ' ... Done. Duration = ' + str(time.clock() - starting_time) + '\n')
    sys.stdout.flush()

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
        description='Caches the top cells.'
    )
    
    # Provide all arguments.
    # Note: we provide nargs=N and the N arguments from the command line will be gathered together into a list.
    # Thus, we supress nargs.
    parser.add_argument('-hh',           required=True,  action='store', type=int, dest='h',             metavar='arg',  help='All top cells with h = 2g+m-1 = 2*genus + number of outgoing boundaries - 1.')
    parser.add_argument('-p', '--parametrized',          action='store_true',      dest='parametrized',                  help='Consider parametrized boundary.', default=False)
    args=vars( parser.parse_args() )
    
    # Create directories.
    pyradbar.create_directories(0, 0)

    # Write Preamble.
    pre, valid = pyradbar.preamble()
    sys.stdout.write(pre + '\n')
    sys.stdout.flush()
    
    # Stop of something went wrong.
    if valid == False:
        print "Could not initialize everything. Abroating." 
        return 1
    
    # Start actual computation.
    for h in range( 1,  args['h']+1 ):
        store_top_cells( h,  args['parametrized'] )
    
    # Cleanup screen with new lines. This is not written to the result file.
    sys.stdout.write('\n\n\n')
    sys.stdout.flush()

if __name__ == "__main__":
    main()
