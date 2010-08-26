#!/usr/bin/env python
#
# pKaTool - analysis of systems of titratable groups
# Copyright (C) 2010 Jens Erik Nielsen
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact information: 
# Email: Jens.Nielsen_at_gmail.com
# Normal mail:
# Jens Nielsen
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland

A={}
A[0]={0:0,
      1:1.3,
      2:2.25,
      3:3.11,
      4:3.42}

A[1]={1:0,
      2:1.3,
      3:2.25,
      4:2.85}

A[2]={2:0,
      3:1.3,
      4:2.25}

A[3]={3:0,
      4:1.3}

A[4]={4:0}

for x in A.keys():
    for y in A[x].keys():
        if not A[y].has_key(x):
            A[y][x]=A[x][y]
x=A.keys()
x.sort()
print x
import dist_geom
X=dist_geom.distance_optimisation(A)
