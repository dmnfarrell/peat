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

import sys,os, pickle

fd=open(sys.argv[1])
dict=pickle.load(fd)
fd.close()

print 'Dict loaded'
mkeys=dict.keys()
mkeys.sort()

#
# Test we have all
#
if dict.has_key('pdbfile'):
    print 'PDB file present'
else:
    print 'PDB file missing'

#
# Get the groups
#
if not dict.has_key('wt_full'):
    print 'wild type pKa values are missing'
    raise 'pKa values missing'

groups=dict['wt_full'].keys()
groups.sort()
print '# of titratable groups: %d' %(len(groups))

for group in groups:
    text='Completeness of %15s ...' %group
    incomplete=None
    print text,
    sys.stdout.flush()
    if dict.has_key(group):
        sub=dict[group]
        cont=None
        if not (sub.has_key('+2.0') and sub.has_key('m2.0')):
            kkeys=sub.keys()
            newkeys=[]
            for key in kkeys:
                if key!='locked':
                    newkeys.append(key)
            text= 'Incomplete design phase:%s' %(str(newkeys))
            print text
            for k2 in kkeys:
                if k2=='locked':
                    continue
                print k2,sub[k2].keys()
            #print kkeys,
            cont=1
        if sub.has_key('locked'):
            if sub['locked']==1:
                print 'Group is LOCKED',
                cont=1
            else:
                if cont==1:
                    if not sub.has_key('pKa out of range'):
                        print 'Group is unlocked and incomplete ???',
                        cont=1
                
        if cont:
            print
            continue
        #
        # Examine the solutions
        #
        for crit in sub.keys():
            if crit=='locked':
                continue
            sol=sub[crit]
            if sol:
                solutions=sol.keys()
                solutions.sort()
                text='%6s:%3d solutions' %(crit,len(solutions))
                print text,
                for solution in solutions:
                    if len(sol[solution].keys())!=6:
                        print sol[solution].keys()
                        incomplete=1
            else:
                incomplete=1
                text='%6s:None' %crit
                print text,
    else:
        incomplete=1
    if incomplete:
        print '.........Incomplete!!!'
    else:
        print '.........Complete'
