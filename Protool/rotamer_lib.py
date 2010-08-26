#!/usr/bin/env python
#
# # Protool - Python class for manipulating protein structures
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


"""
Read the rotamer library
"""
import Protool, os, cPickle
dir_location=os.path.split(Protool.__file__)[0]
location=os.path.join(dir_location,'bbdep02.May.sortlib')
text='Importing rotamer library for Protool from %s.....' %location
print text,
import sys, mmap

sys.stdout.flush()
local_pickle=os.path.join(dir_location,'bbdep02.pickle')

#
# Parse the original file
#
fd=open(location)
line=fd.readline()
rots={}
while line:
	sp=line.split()
	name=sp[0]
	if not rots.has_key(name):
	    rots[name]=[]
	#
	# Fill in the data
	#
	tr={}
	tr['phi']=float(sp[1])
	tr['psi']=float(sp[2])
	tr['p']=float(sp[8])
	tr['chi1']=float(sp[9])
	tr['chi2']=float(sp[10])
	tr['chi3']=float(sp[11])
	tr['chi4']=float(sp[12])
	#
	# Add this rotamer
	#
	rots[name].append(tr.copy())
	#
	# Next line
	#
	line=fd.readline()

fd.close()
#
# Rotamer library read
#
print 'Done'

