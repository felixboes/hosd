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

import itertools

class Permutation(object):

  def __init__(self, d, perm = None):
    self._d = d
    self._norm = -1
    if perm is None:
      self._rep = tuple(range(d + 1))
    else:
      self._rep = tuple(perm)
    self._hash = None

  def compute_norm(self):
    i = 0
    self._norm = self._d
    visited = (self._d + 1) * [False]
    visited[self._d] = True

    while not visited[i]:
      self._norm -= 1
      j = self._rep[i]
      visited[j] = True
      while j != i:
        j = self._rep[j]
        visited[j] = True

      while visited[i] and i < self._d:
        i += 1

  def degree(self):
    return self._d - 1

  def norm(self):
    if self._norm == -1: # Norm has to be computed.
      self.compute_norm()
    return self._norm

  def num_cyc(self):
    return self._d - self.norm()

  def fixed_pts(self):
    return [x for x in range(self._d) if self._rep[x] == x]

  def num_fixed_pts(self):
    return len(self.fixed_pts())

  def cycle_decomposition(self):
    i = 0
    cycle_decomp = []
    visited = (self._d + 1) * [False]
    visited[self._d] = True

    while not visited[i]:
      cycle_decomp.append([i])
      j = self._rep[i]
      visited[j] = True
      while j != i:
        cycle_decomp[-1].append(j)
        j = self._rep[j]
        visited[j] = True

      while visited[i] and i < self._d:
        i += 1

    return tuple(tuple(x) for x in cycle_decomp)

  def face(self, num_i):
    i = int(num_i) % self._d

    suc = self._rep[i]
    transp = [x for x in range(self._d)]
    transp[i] = suc
    transp[suc] = i
    intermediate_tau = [transp[self._rep[x]] for x in range(self._d) if x != i]

    return Permutation(self._d - 1, [x if x < i else x - 1 for x in intermediate_tau] )

  def __str__(self):
    return str(self.cycle_decomposition())

  def __eq__(self, other):
    return isinstance(other, Permutation) and self._d == other._d and self._rep == other._rep

  def __hash__(self):
    if self._hash is None:
      self._hash = int(0)
      for x in range(self._d):
        self._hash += self._rep[x] * 10 ** x
      return self._hash
    else:
      hash = int(0)
      for x in range(self._d):
        hash += self._rep[x] * 10 ** x
      if self._hash == hash:
        return self._hash
      else:
        raise RuntimeError('Hash values are different...', self._hash, ' ', hash, sep='')


class SymmetricGroup(object):

  def __init__(self, d):
    self._d = d

  def get_long_cycles(self):
    # Compute permtuations that have one cycle by giving the cycle decomposition:
    cyc_decs = [(0,) + p for p in itertools.permutations( range(1, self._d) )]
    # Compute the map representation:
    permutations = []
    for cyc in cyc_decs:
      rep = self._d*[0]
      for i in range(self._d):
        rep[cyc[i]] = cyc[(i + 1) % self._d]
      permutations.append(Permutation(self._d, rep))
    return permutations