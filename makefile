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

LANG = en_US.UTF-8

CXXFLAGS      := -O3 -std=c++11 -Wextra -Wall -Wno-long-long -pedantic-errors
INCL          := -I./chomp/include
LIBS          := -lgmpxx -lgmp
SRC_DIR       := ./chomp

ifndef CXX
CXX           := g++
endif

.PHONY: all
all: homchain homchain_gmp

.PHONY: homchain_gmp
homchain_gmp:
	$(CXX) $(CXXFLAGS) -DCHOMP_GMP_VERSION=1 $(INCL) -o homchain_gmp $(SRC_DIR)/programs/homprogs/homchain.cpp $(SRC_DIR)/src/system/arg.cpp $(SRC_DIR)/src/system/textfile.cpp $(SRC_DIR)/src/system/timeused.cpp $(SRC_DIR)/src/struct/integer.cpp $(LIBS)

.PHONY: homchain
homchain:
	$(CXX) $(CXXFLAGS) $(INCL) -o homchain $(SRC_DIR)/programs/homprogs/homchain.cpp $(SRC_DIR)/src/system/arg.cpp $(SRC_DIR)/src/system/textfile.cpp $(SRC_DIR)/src/system/timeused.cpp $(SRC_DIR)/src/struct/integer.cpp $(LIBS)

.PHONY: clean
clean:
	rm -rf homchain homchain_gmp
