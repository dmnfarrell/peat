#!/usr/bin/env python
#
# pKarun - scripts for running pKa calculations
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

import os

def rmr(dir):
    """Remove a dir recursively (rm -r) or if it is a file: unlink it"""
    if os.path.isdir(dir):
        for file in os.listdir(dir):
            if os.path.isdir(os.path.join(dir,file)):
                rmr(os.path.join(dir,file))
            if os.path.isfile(os.path.join(dir,file)) or os.path.islink(os.path.join(dir,file)):
                os.unlink(os.path.join(dir,file))
        if os.path.isdir(dir):
            #print 'RMDIR:',dir
            os.rmdir(dir)
        return
    else:
        os.unlink(dir)
        return

def delete(list):
    """Delete a file if it is there"""
    for file in list:
        if os.path.isfile(file):
            os.unlink(file)
    return

def find(pattern, dir=os.getcwd()):
    import fnmatch
    list = []
    names = os.listdir(dir)
    names.sort()
    for name in names:
            fullname=os.path.join(dir,name)
            if os.path.isfile(fullname):
                if fnmatch.fnmatchcase(name,pattern):
                    list.append(fullname)
            elif os.path.isdir(fullname):
                if fnmatch.fnmatchcase(name,pattern):
                    list.append(fullname)
                #print 'Descend into ',fullname
                list=list+find(pattern,fullname)
    return list                                               
