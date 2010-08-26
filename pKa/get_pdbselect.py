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

pdbselect='/enzyme/nielsen/pKD_bench/recent.pdb_select25'
destination='/home/nielsen/pKa-design/pdbselect'

import os
import string

files={}
fd=open(pdbselect)
lines=fd.readlines()
fd.close()
failed=[]
total=len(lines)
count=0
for line in lines:
    #
    # Keep a tab on how far we've gotten
    #
    pdbid=string.split(line)[1]
    print pdbid,total-count
    count=count+1
    #
    # Start the work
    #
    id=string.lower(pdbid[:4])
    keep_chainid=None
    #print 'id',id
    if len(pdbid)>4:
        keep_chainid=pdbid[4]
        if keep_chainid=='_':
            keep_chainid=None
    #
    # Set all filenames
    #
    dir=os.path.join(destination,pdbid)
    pdb=os.path.join('/data/pdb',id+'.pdb')
    newpdb=os.path.join(dir,id+'.pdb')
    pkapdb=os.path.join(dir,id+'.pka.pdb')
    if os.path.isfile(pkapdb):
        continue
    #
    # Make new dir
    #
    if not os.path.isdir(dir):
        os.mkdir(dir)
    #
    # Does the PDB file exist?
    #
    if not os.path.isfile(pdb):
        failed.append(pdbid)
        continue
    os.system('cp %s %s' %(pdb,newpdb))
    #
    # Correct with WI
    #
    import WI_tools
    renumb=1
    setcha=None
    readallmodels=None
    logifile,files=WI_tools.delwat_dellig_corall(newpdb,renumb,readallmodels,setcha)
    filenms=files.keys()
    if len(filenms)!=1:
        print 'More than one filename',filenms
        failed.append(pdbid)
        continue
    import string
    #
    # Do we have a chain identifier?
    #
    chainid=-1
    changes=0
    for line in files[filenms[0]]:
        split=string.split(line)
        #print split
        if split[0]=='ATOM':
            if line[21]!='':
                chainid=1
            #if split[4] in string.letters and len(split[4])==1:
            #    if not chainid:
            #        changes=changes+1
            #        chainid=1
            #else:
            #    if chainid==1:
            #        changes=changes+1
            #        chainid=None
            #    elif chainid==-1:
            #        chainid=None
    #if changes>1:
    #    raise 'Could not determine if we have a chainid or not',pdbid
    #
    # Write the PDB file
    #
    #print 'chainid',chainid
    #print 'keepchainid',keep_chainid
    fd=open(pkapdb,'w')
    for line in files[filenms[0]]:
        split=string.split(line)
        newline=line
        if split[0]=='ATOM' and keep_chainid:
            if chainid:
                thisid=split[4]
                thisid=line[21]
                #print string.upper(thisid),string.upper(keep_chainid)
                if string.upper(thisid)!=string.upper(keep_chainid):
                    newline=''
                    #print 'Excluding line'
                #else:
                #    print 'Keeping',line
        elif split[0]=='HETATM':
            newline=''
        fd.write(newline)
    fd.close()
    #
    # Remove the chain ID
    #
    import remove_chain_identifier
    remove_chain_identifier.remove_chainid(pkapdb)
    #
    # Done
    #
print
print 'Failed:'
print failed
