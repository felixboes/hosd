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

def norm( perm ):
    p = len( perm )
    i = 0
    norm = p
    visited = (p+1)*[False]
    visited[p] = True
    
    while( visited[i] == False ):
        norm -= 1
        j = perm[i]
        visited[j] = True
        while( j != i ):
            j = perm[j]
            visited[j] = True
        
        while( visited[i] == True and i < p ):
            i += 1
    
    return norm

def num_cyc( perm ):
    return len(perm) - norm(perm)
    
def fixed_pts( perm ):
    p = len( perm )
    return [ x for x in range(p) if perm[x] == x ]

def num_fixed_pts( perm ):
    return len( fixed_pts(perm) )

def cycle_decomposition( perm ):
    p = len( perm )
    i = 0
    cycle_decomp = []
    visited = (p+1)*[False]
    visited[p] = True
    
    while( visited[i] == False ):
        cycle_decomp.append( [i] )
        j = perm[i]
        visited[j] = True
        while( j != i ):
            cycle_decomp[-1].append(j)
            j = perm[j]
            visited[j] = True
        
        while( visited[i] == True and i < p ):
            i += 1
    
    return tuple( tuple(x) for x in cycle_decomp )

def d(tau, num_i):
    p = len(tau)
    i = int(num_i) % p
    
    suc = tau[i]
    transp = [ x for x in range(p) ]
    transp[ i ] = suc
    transp[ suc ] = i
    intermediate_tau = [ transp[ tau[x] ] for x in range(p) if x != i  ]
    
    return tuple ( x if x < i else x - 1 for x in intermediate_tau )
