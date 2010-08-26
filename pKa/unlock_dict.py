#!/usr/bin/env python
#
# pKa - various programs and scripts for pKa value analysis, calculation and redesign
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


#
# Unlocks the dictionary for Design_accuracy after a crash
#

import sys
fd=open(sys.argv[1])
import pickle
dict=pickle.load(fd)
fd.close()
import types
for key in dict.keys():
    if type(dict[key]) is types.DictType:
        if dict[key].has_key('locked'):
            print 'Unlocking',key
            del dict[key]['locked']

fd=open(sys.argv[1],'w')
pickle.dump(dict,fd)
fd.close()

import os
lockfile=sys.argv[1]+'.lock'
if os.path.isfile(lockfile):
    os.unlink(lockfile)
