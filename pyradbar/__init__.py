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

import sys
major, minor = sys.version_info[0], sys.version_info[1]
if major == 2:
  from cell import *
  from cell_para import *
  from draw_cell import *
  from misc_useability_stuff import *
  from perm import *
  from symmetric_group import *
elif major == 3:
  from .cell import *
  from .cell_para import *
  from .draw_cell import *
  from .misc_useability_stuff import *
  from .perm import *
  from .symmetric_group import *
  from .homology_computation_symm import compute_homology_symm
