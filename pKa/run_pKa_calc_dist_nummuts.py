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

import os, string
pwd=os.getcwd()
files=os.listdir(pwd)
pdbfile=None
for file in files:
    if string.find(file,'pka.pdb')!=-1 and string.find(file,'pka.pdb.')==-1:
        pdbfile=file
        break
if not pdbfile:
    print 'Could not find pdb file???'
#
# Check that we didn't do this already
#
pkafile=pdbfile+'.PKA.DAT'
#
# Got the PDB file
#
# Make symlinks
#
source='/enzyme/nielsen/pkaparms'
files=['DELRAD.DAT','DELCRG.DAT','TOPOLOGY.H']
for file in files:
    r_s=os.path.join(source,file)
    r_d=os.path.join(pwd,file)
    os.system('ln -s %s %s' %(r_s,r_d))
#
# Start pKa calc
#
import pKa
params={'dbcrit':1000,'lowph':0.1,'highph':20.0,'phstep':0.1,'pairene':10.0}
Y=pKa.pKarun(os.getcwd(),pdbfile,params)
if not os.path.isfile(pkafile):
    Y.runseq()
#
# Do the design
#
import Design_dist_nummuts
Design_dist_nummuts.main(pdbfile)
#
# All done
#
