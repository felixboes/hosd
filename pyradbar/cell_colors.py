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

def cell_color( idx ):
    cell_colors = { \
       -1 : [1, 1, 1], \
        0 : [0, 0, 0], \
        1 : [1, 0, 0], \
        2 : [0, 1, 0], \
        3 : [0, 0, 1] \
    }
    
    cell_colors_by_name = { \
        'white' : [1, 1, 1], \
        'black' : [0, 0, 0], \
        'red'   : [1, 0, 0], \
        'green' : [0, 1, 0], \
        'blue'  : [0, 0, 1]
    }
    
    if isinstance( idx, int ) :
        return cell_colors[ idx % ( len(cell_colors) - 1 ) ]
    elif isinstance( idx,  basestring):
        return cell_colors_by_name[ idx ]
