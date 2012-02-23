#!/usr/bin/env python
#
# DataPipeline - A data import and fitting tool
# Copyright (C) 2011 Damien Farrell
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
# Email: damien.farrell_at_ucd.ie
# Normal mail:
# Damien Farrell
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland
#

'''Module containing utility classes and functions'''

import os

def setAttributesfromConfigParser(obj, cp):
    """A helper method that makes the options in a ConfigParser object
       attributes of obj"""
   
    for s in cp.sections():
        obj.__dict__[s] = cp.items(s)          
        #print cp.items(s)
        for f in cp.items(s):
            #print f[0], f[1]
            try: val=int(f[1])
            except: val=f[1]
            obj.__dict__[f[0]] = val

def getListFromConfigItems(items):
    """Get a list from a set of ConfigParser key-value pairs"""
    lst = [i[1] for i in sorted(items) if i[1] != '']        
    return lst

def clearDirectory(path):
    """Remove all files in folder"""
    for f in os.listdir(path):
        filepath = os.path.join(path, f)
        try:
            if os.path.isfile(filepath):
                os.unlink(filepath)
        except Exception, e:
            print e
    return    
