#!/usr/bin/env python
#
# pKaTool - analysis of systems of titratable groups
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

def main():
    import sys
    file1=sys.argv[1]
    file2=sys.argv[2]
    import pKaTool
    X=pKaTool.pKaIO.pKaIO()
    pka1=X.readpka(file1)
    pka2=X.readpka(file2)
    #
    # Loop over all residues and calculate rmsd
    #
    residues1=pka1.keys()
    residues2=pka2.keys()
    residues1.sort()
    residues2.sort()
    if residues1!=residues2:
        print 'The two files do not contain the same titratable groups'
        raise Exception()
    diffs=[]
    for res in residues1:
        if pka1[res]['pKa'] and pka2[res]['pKa']:
            print '%s pka1: %5.2f, pka2 %5.2f, diff: %5.3f' %(res,pka1[res]['pKa'],pka2[res]['pKa'],abs(pka1[res]['pKa']-pka2[res]['pKa']))
            diffs.append(abs(pka1[res]['pKa']-pka2[res]['pKa']))
        else:
            print 'Skipping %s. pka1 %s, pka2 %s' %(res,pka1[res]['pKa'],pka2[res]['pKa'])
    sumsq=0.0
    for diff in diffs:
        sumsq=sumsq+diff*diff
    sumsq=sumsq/float(len(diffs))
    import math
    print 'RMSD: %5.3f' %(math.sqrt(sumsq))
    return

if __name__=='__main__':
    main()
