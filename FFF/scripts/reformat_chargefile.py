#!/usr/bin/env python
#
# FFF - Flexible Force Field
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
import sys
crgfile=sys.argv[1]
radfile=sys.argv[2]

def get_values(filename):
    fd=open(filename)
    lines=fd.readlines()
    fd.close()

    data={}
    count=0
    record=False
    while 1:
        line=lines[count].strip()
        if len(line)==0:
            count=count+1
            continue
        if line[0]=='*':
            record=True
        if lines[count].strip()=='#END':
            break
        count=count+1
        if not record:
            continue
        #
        done=False
        recordname=lines[count].strip()
        count=count+1
        atoms={}
        #print 'Got record name',recordname
        while not done:
            if len(lines[count].strip())==0:
                done=True
                break
            atom=lines[count].strip().split()
            atoms[atom[0]]=float(atom[1])
            count=count+1
        data[recordname]=atoms.copy()
    return data
    
if __name__=='__main__':
    crg=get_values(crgfile)
    rad=get_values(radfile)
    print '#BEGIN'
    records=crg.keys()
    records.sort()
    for record in records:
        print '*'
        print record
        atoms=crg[record].keys()
        atoms.sort()
        for atom in atoms:
            charge=crg[record][atom]
            radius=rad[record][atom]
            print '%5s %6.3f %6.3f' %(atom,charge,radius)
        print
    print '#END'
            
    