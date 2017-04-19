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

import copy
import os
import sys
import gzip
if sys.version_info[0] == 2:
    import cPickle as pickle
    import perm
else:
    import pickle
    from . import perm
import inspect


class Cell:
    def __init__( self, rho, suspend=True ):
        # Right now, we can create top cells with g,m != 0,1 correctly.
        # This is ok, since we do not need more.
        self._p = len(rho)+1 if suspend else len(rho)
        if self._p > 1:
            p_2 = self._p-2
        else:
            p_2 = 0
        if suspend:
            # rho[x] = p-2 - > rho[x] = p-1
            # rho[p-1] = p-2
            self._rho = tuple( [ (rho[x] if rho[x] != p_2 else self._p-1) for x in range(self._p-1)] ) + tuple( [p_2] )
            # shift by one.
            self._rho = tuple( [ (self._rho[(x+self._p -1)%self._p]+1)%self._p for x in range(self._p) ] )
        else:
            self._rho = rho
        self._lam = None
        self._lam_cyc_of_symbol = None
        self.compute_lam()
        self._m = self.num_cyc_rho() + self.num_fixed_pts_lam() - (1 if suspend else 0)
        # An Euler characteristic argument shows (for desuspended cells)
        # 2 - 2g - m - 1 = - num pixed pts lam - norm lam; so g = (num pixed pts lam + norm lam + 1 - m)/2
        self._g = ( self.num_fixed_pts_lam() - (1 if suspend else 0) + self.norm_lam() - self._m + 1 ) /2
        cycles_of_lam = self.cyc_lam()
        self._surface_of_symbol = copy.deepcopy(self._lam_cyc_of_symbol)
        self._surface = [ [ [n], 0, 0 ] for n, cyc in enumerate(cycles_of_lam) ]
        self._degenerated = False
        self._hash = None
        
        self._valid = True
    
    def compute_lam(self):
        # lam is (0,1,2,...,p) * rho, i.e.
        # x maps to rho(x) + 1 mod p
        self._lam = tuple( (self._rho[x] + 1) % self._p for x in range(self._p) )
        self._lam_cyc_of_symbol = self._p*[0]
        i = 0
        for cyc in self.cyc_lam():
            for sym in cyc :
                self._lam_cyc_of_symbol[sym] = i
            i += 1
    
    def norm_rho(self):
        return perm.norm( self._rho )
    
    def cyc_rho(self):
        return perm.cycle_decomposition( self._rho )
    
    def num_cyc_rho(self):
        return perm.num_cyc( self._rho )
    
    def norm_lam(self):
        return perm.norm( self._lam )
    
    def cyc_lam(self):
        return perm.cycle_decomposition(self._lam)
    
    def num_cyc_lam(self):
        return perm.num_cyc( self._lam )
    
    def fixed_pts_lam(self):
        return perm.fixed_pts( self._lam )
    
    def num_fixed_pts_lam(self):
        return perm.num_fixed_pts( self._lam )
    
    def num_surf(self):
        return len(self._surface)
    
    def degree(self):
        return self._p - 1
    
    def get_rho(self):
        return self._rho
    
    def get_lam(self):
        return self._lam
    
    def get_num_punc(self):
        return self._m
    
    def get_genus(self):
        return self._g
    
    def get_degenerated(self):
        return self._degenerated
    
    def get_surfaces(self):
        return self._surface
    
    # In Python 3.* we have a __bool__() function.
    # In Python 2.* we are forves to use __nonzero__().
    def __bool__(self):
        return self._valid
    __nonzero__ = __bool__
    
    def face(self, num_i):
        # check degenerateness
        if self._p == 1:
            self._valid = False
            return
        
        # copy old status
        old_p = copy.deepcopy(self._p)
        old_i = num_i % old_p
        old_lam = copy.deepcopy( self._lam )
        old_lam_cyc_of_symbol = copy.deepcopy( self._lam_cyc_of_symbol )

        # compute new status but be warned,
        # self._surface and self._surface_of_symbol are computed during the process.
        difference = self.num_cyc_lam()
        self._rho = perm.d(self._rho, num_i)
        self._p -= 1
        self.compute_lam()
        i = num_i % self._p
        difference = self.num_cyc_lam() - difference
        
        # case analysis how num_cyc_lam changes
        # -1 iff two boundary circles of surfaces collide
        # +0 iff puncture degenerates
        # +1 iff a boundary circle is split into two
        if difference == -1 :
            cyc_idx_1 = old_lam_cyc_of_symbol[ old_i ]
            sur_idx_1 = self._surface_of_symbol[ old_i ]
            cyc_idx_2 = old_lam_cyc_of_symbol[ (old_i + 1) % (old_p) ]
            sur_idx_2 = self._surface_of_symbol[ (old_i + 1) % (old_p) ]
            # Since cycles are ordered wrt their smalles symbol in their support,
            # the index of the new cycle is the minimum.
            # We overwrite the old cycle at this position and remove the cycle at the maximum.
            cyc_idx_new = min(cyc_idx_1, cyc_idx_2)
            cyc_idx_rem = max(cyc_idx_1, cyc_idx_2)
            
            # We discuss case of two colliding boundary cycles of the same surface first.
            # In any case, we have to renormalize the cycle numbers.
            if sur_idx_1 == sur_idx_2:
                # Here, the genus is increased by one. 
                self._surface[ sur_idx_1 ][1] += 1
                self._degenerated = True
            else:
                # Here, two surfaces are connected.
                # As above, surfaces are ordered wrt their cycles
                sur_idx_new = min(sur_idx_1, sur_idx_2)
                sur_idx_rem = max(sur_idx_1, sur_idx_2)
                
                # The new surface is the sum of the old two.
                # Since we renormalize later, we may keep the cycle that is going to be removed for a while.
                surf_new = [ self._surface[ sur_idx_new ][j] + self._surface[ sur_idx_rem ][j] for j in range(3) ]
                surf_new[0].sort()
                
                # We overwrite the one of the old surfaces with the new one and remove the other.
                self._surface[sur_idx_new] = surf_new
                self._surface.pop(sur_idx_rem)
                
                # We have to recompute the map which associates to a symbol a surface.
                self._surface_of_symbol = [ x if x != sur_idx_rem else sur_idx_new for x in self._surface_of_symbol ]
                self._surface_of_symbol = [ x if x < sur_idx_rem else x - 1 for x in self._surface_of_symbol ]
            
            # The cycle with index cyc_idx_rem is removed, the others are renormalized.
            for j, surf in enumerate(self._surface):
                self._surface[j] = [ [ x if x < cyc_idx_rem else x - 1 for x in surf[0] if x != cyc_idx_rem ], surf[1], surf[2] ]
        
        elif difference == 0:
            # increase weight by one.
            self._surface[ self._surface_of_symbol[ i ] ][2] += 1
            self._degenerated = True
        else:
            # The surface splits a boundary cycle into two parts.
            # A picture will show clearify the role of a.
            old_a = old_lam[old_i]
            a = old_a - 1 if old_i < old_a else old_a
            
            # The new indices are given by i and a.
            cyc_idx_1 = self._lam_cyc_of_symbol[ i ]
            cyc_idx_2 = self._lam_cyc_of_symbol[ a ]
            
            # The cycle that stays in its position is cyc_idx_old.
            cyc_idx_old = min(cyc_idx_1, cyc_idx_2)
            cyc_idx_new = max(cyc_idx_1, cyc_idx_2)
    
            # Renormalize the cycles.
            for j, surf in enumerate(self._surface):
                self._surface[j] = [ [ x if x < cyc_idx_new else x + 1 for x in surf[0] ], surf[1], surf[2] ]
            # Put the new cycle into the surface.
            self._surface[ self._surface_of_symbol[i] ][0].append( cyc_idx_new )
            self._surface[ self._surface_of_symbol[i] ][0].sort()
        
        # Remove the symbol i from the map associating a surface to a symbol.
        self._surface_of_symbol = [ x for j, x in enumerate(self._surface_of_symbol) if j != old_i ]
    
    def orientation_sign(self, num_i):
        # Compute to which cycle each symbol belongs.
        num_cycles = 0
        cycles_of_symbol = self._p*[0]
        cycles_of_rho = self.cyc_rho()
        for cycle in cycles_of_rho:
            for symbol in cycle :
                cycles_of_symbol[symbol] = num_cycles
            num_cycles += 1
        
        cycle_num_of_i = cycles_of_symbol[num_i]
        cycle_of_i = cycles_of_rho[ cycle_num_of_i ]
        
        # If i is not the smallest element in its cycle, then the order of the cycles remain unchanged.
        # Else we remove i (which was the smallest symb) and compute the new position of the cycle.
        delta = 0
        if num_i == min( cycle_of_i ):
            if len(cycle_of_i) == 1:
                return 0
            j = min ( x for x in cycle_of_i if x != num_i )
            delta = 0
            for k in range( num_i + 1, j ):
                if cycle_num_of_i + delta < cycles_of_symbol[k]:
                    delta += 1
        return -1 if delta % 2 == 0 else 1
    
    def rotation_symmetric(self):
        # A cell is considered rotation symmetric if the resp. desuspended
        # cell is.
        cell = self.get_clean_copy()
        if cell._lam[0] == 0:
            cell.face(0)
        
        # A cell is rotation symmetric if and only if 
        # lambda is an affine function from Z_p to Z_p.
        offset = cell._lam[0]
        for i in range(0, cell._p):
            if cell._lam[i] != (i + offset) % cell._p:
                return False
        
        return True
    
    def starts_cycle(self):
        # We assert:
        #    The cell is essential wrt first Morseflow.
        #
        # The function does the following.
        # 1) We make the cell desuspended and make all odd fans even.
        #    This is reasonable, since a cell produces in a cycle
        #    if and only if this version of the cell does.
        # 2) We rotate the cell (by shifting all symbols by one)
        #    and extract the minimal word out of the correponding words
        #    (lambda(0), .... lambda(p-1))
        #    The cell is considered in minimal position of the minimal word is
        #    given by the cell.
        
        cell = self.get_clean_copy()
        # A suspended cell is in minimal position if and only if
        # its desuspended face is.
        if cell._lam[0] == 0:
            cell.face(0)
        
        # Now, we make all odd cells even.
        i = 0
        do_while = True
        while do_while:
            fan_start = i
            
            # Compute the length of the current fan.
            while cell._lam[i] == i+1:
                i = i+1
            
            # Check if the do_while loop has to end.
            if i + 1 == cell._p:
                do_while = False
                # Check if the fan ends at zero or p-1
                if cell._lam[i] == 0:
                    i = i+1
            
            # Check if the fan is odd.
            if (i - fan_start) % 2 == 1:
                # If the odd fan is isolated, the cell is not fully rotatable hence it
                # cannot start a cycle.
                # Here, we use that the cell is essential wrt first Morse flow.
                # Dropping this assertion, we do not properly detect odd fans at zero.
                if i != cell._p:
                    # This is an isolated odd fan.
                    if cell._lam[ i ] == fan_start:
                        return False
                    else:
                        # Otherwise, we make the fan even.
                        cell.face(fan_start)
                        i = i-1
                else:
                    # This is an isolated odd fan under the condition,
                    # that the cell is not rotation symmetric.
                    if cell._lam[ 0 ] == fan_start or cell._lam[ 0 ] == 1:
                        return False
                    else:
                        # Otherwise, we make the fan even.
                        cell.face(fan_start)
                        i = i-1
            
            # Next i
            i = i+1
        
        # Iterate through all conjugates
        # (0 ... k)^l lam (0 ... k)^{-l}
        # and take the minimal word.
        minimal = copy.deepcopy(cell._lam)
        for l in range(1, cell._p):
            rotated = tuple( [ (cell._lam[ (x - l + cell._p) % cell._p ] + l) % cell._p for x in range(cell._p) ] )
            if rotated < minimal:
                minimal = rotated
        
        if cell._lam == minimal:
            return True
        else:
            return False
    
    def morse_type(self):
        # The single cell in degree zero is essential.
        if self.degree() == 0:
            return 'e', None
        
        # In order to iterate through the minimal fans,
        # we iterate through the cycles of lambda.
        cycles_of_lam = self.cyc_lam()
        
        for number_cyc, cyc in enumerate(cycles_of_lam):
            # This is the minimal fan of the current cycle if g = 0.
            # For g > 0, we have to check the minimality condition.
            fan_start = cyc[0]
            if number_cyc == min( self._surface[ self._surface_of_symbol[fan_start] ][0] ):
                # Compute its length.
                i = fan_start
                while self._lam[ i ] == i+1 :
                    i = i+1
                len_fan = i - fan_start
                
                # The case of the long fan has to be handled differently.
                if fan_start == 0 and i+1 == self._p:
                    len_fan += 1
                
                # The cell is collapsible if the length is odd.
                if len_fan % 2 == 1:
                    return 'c', fan_start
                # The cell is redundant if the length is even and
                # there are extendible punctures.
                elif self._surface[self._surface_of_symbol[fan_start]][2] > 0:
                    return 'r', None
        
        # Now, we are going to enlarge the Morse flow.
        cell_is_suspended = True if self._lam[0] == 0 else False
        
        # A rotation symmetric cell has cannot start an oriented cycle.
        # One problem is the invertablity of the coefficient:
        # It is either zero or two.
        if not self.rotation_symmetric() and not self.starts_cycle():
            if cell_is_suspended:
                return 'c', 0
            else:
                return 'r', None
        
        return 'e', None
    
    def get_filename_str(self):
        return ( 'g=' + str(self._g) + ',m=' + str(self.get_num_punc()) + ',lam=' + str(self.cyc_lam()) ).replace(" ", "")
    
    def get_clean_copy(self):
        ret = copy.deepcopy(self)
        ret._hash = None
        return ret
    
    def __str__(self):
        if( self._valid == True ):
            return ''.join( str(e) for e in [
                "Genus = ", self._g, "; ", 
                "Number of Punctures = ", self._m, "; ", 
                "Rho = ", self.cyc_rho(),  "; ", 
                "Lambda = ", self.cyc_lam(),  "; ", 
                "Surfaces = " , self._surface ] )
        else:
            return 'Invalid Cell: ' + str( self._rho )
    
    def __hash__(self):
        if self._hash is None:
            self._hash = int(0)
            for x in range(len(self._rho)):
                self._hash += self._rho[x]*10**x
            return self._hash
        else:
            hash = int(0)
            for x in range(len(self._rho)):
                hash += self._rho[x]*10**x
            if self._hash == hash:
                return self._hash
            else:
                sys.stdout.write("Hash values are different..." + str(self._hash) + " " + str(hash) + "\n")
                return 0
    
    def __eq__(self, other):
        if self._rho == other._rho and self._surface == other._surface :
            return True
        else:
            return False

class TopCellRhoArchive:
    def __init__(self, gziparchive=None):
        self._archive = gziparchive
    
    def __iter__(self):
        return self
    
    def next(self):
        if self._archive is None:
            raise StopIteration
        
        try:
            rho = pickle.load( self._archive )
        except EOFError:
            # We reached the archives end.
            raise StopIteration
        except:
            # Other error.
            frameinfo = inspect.getframeinfo(inspect.currentframe())
            sys.stdout.write(str(frameinfo.filename) + ' ' + str(frameinfo.lineno) + '\n')
            e, p, t = sys.exc_info()
            sys.stdout.write( str(e) + ' ' + str(p) + '\n')
            raise StopIteration
        else:
            return rho
    # Python 3:
    __next__ = next
    
    def close(self):
        if isinstance(self._archive, gzip.GzipFile):
            self._archive.close()
        elif self._archive is not None:
            sys.stdout.write("Wrong file type: " + str(type(self._archive)) + '\n')

class LoadTopCellRho:
    def __init__(self, g, m):
        self._g = g
        self._m = m
        self._valid = True
        self._archive = None
        if g < 0 or m < 1:
            self._valid = False
        
    def __enter__(self):
        # open archive
        try:
            if( os.path.exists('./data/') and os.path.exists('./data/top_cell_g_' + str(self._g) + '_m_' + str(self._m) + '_unparametrized.bz2') ):
                self._archive = TopCellRhoArchive( gzip.GzipFile('./data/top_cell_g_' + str(self._g) + '_m_' + str(self._m) + '_unparametrized.bz2', 'rb') )
                return self._archive
            else:
                sys.stdout.write('File ' + './data/top_cell_g_' + str(self._g) + '_m_' + str(self._m) + '_unparametrized.bz2 does not name a valid file.\n')
                self._valid = False
                return None
        except:
            # Print the error.
            frameinfo = inspect.getframeinfo(inspect.currentframe())
            sys.stdout.write(str(frameinfo.filename) + ' ' + str(frameinfo.lineno) + '\n')
            e, p, t = sys.exc_info()
            sys.stdout.write(str(e) + ' ' + str(p) + '\n')
            return None
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        
        if exc_type is not None:
            if isinstance(self._archive, TopCellRhoArchive):
                self._archive.close()
            return False
        else:
            if isinstance(self._archive, TopCellRhoArchive):
                self._archive.close()
                return True
            elif self._archive is not None:
                sys.stdout.write("Wrong file type: " + str(type(self._archive)) + '\n')
                return False
            else:
                return True
