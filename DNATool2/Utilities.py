#!/usr/bin/env python
#
# DNATool - A program for DNA sequence manipulation
# Copyright (C) 2012- Damien Farrell & Jens Erik Nielsen
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
# Email: farrell.damien_at_gmail.com 


"""Utility functions for DNATool"""

import os, sys, types, random, string

def createRandomStrings(l,n,chars=None,upper=False):
    """create list of l random strings, each of length n"""

    names = []
    if chars == None:
        chars = string.ascii_lowercase
    #for x in random.sample(alphabet,random.randint(min,max)):
        
    if upper == True:        
        chars = [i.upper() for i in chars]
    for i in range(l):
        val = ''.join(random.choice(chars) for x in range(n))
        names.append(val)
    return names
    
def createDummySequence(l):
    """Create a fake dna seq"""
    seq = createRandomStrings(l,1,'agtc',upper=True)
    return seq
    
def createDummyRestrSites(l=200,n=40):
    """n random rest sites for seq of length l"""
    sites = {}
    names = createRandomStrings(l,6)    
    locs = [random.randint(0, l) for i in range(n)]
    for x,name in zip(locs, names):
        sites[name] = {'position':x}            
    return sites
