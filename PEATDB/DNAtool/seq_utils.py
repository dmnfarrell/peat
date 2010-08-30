#!/usr/bin/env python
#
# Protein Engineering Analysis Tool DataBase (PEATDB)
# Copyright (C) 2010 Damien Farrell & Jens Erik Nielsen
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
# 

import sys, string, sys

#
# ---------
#

def get_base_count (seq, base):
    """Get the number of occurrences of the given base(s) in the sequence."""
    count=0 

    if base=='G':
        for b in seq:
            if b=='G':
                count=count+1
    elif base=='C':
        for b in seq:
            if b=='C':
                count=count+1
    elif base=='T':
        for b in seq:
            if b=='T':
                count=count+1
    elif base=='A':
        for b in seq:
            if b=='A':
                count=count+1
    elif base=='GC':
        for b in seq:
            if b=='G' or b=='C':
                count=count+1    
    elif base=='AT':
        for b in seq:
            if b=='A' or b=='T':
                count=count+1
                
    return count  
    

