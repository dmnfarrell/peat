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

def Nonefix(value):
    if value:
        return '%6.3f' %(value)
    else:
        return 'NA'

def main():
    """Calculate the change in chemical shift due to a full charge on each titratable group"""
    import get_dEF, sys
    method=sys.argv[2] 
    X=get_dEF.map_struct(sys.argv[1])
    X.build_Hs()
    residues=X.PI.residues.keys()
    residues.sort()
    #
    # Get the titratable groups
    #
    titgroups=X.PI.get_titratable_groups()
    titgroups.sort()
    import pKa.pKD_tools
    for titgroup in titgroups:
        titgroup_type=pKa.pKD_tools.get_titgroup_type_from_titgroup(titgroup)
        charge=X.PI.titgroups[titgroup_type]['charge']
        print 'TITRATABLE GROUP',titgroup
        print 'Residue  CS Nitrogen    CS Hydrogen'
        for residue in residues:
            dCS_N=X.get_dCS(residue+':N',titgroup,charge=charge,method=method)
            dCS_H=X.get_dCS(residue+':H',titgroup,charge=charge,method=method)
            print '%8s,    %s,    %s' %(residue,Nonefix(dCS_N),Nonefix(dCS_H))

if __name__=='__main__':
    import sys, os
    error=None
    if len(sys.argv)!=3:
        error=1
    elif not os.path.isfile(sys.argv[1]):
        error=1
    if not error:
        main()
    else:
        print
        print 'Calculate effect of charged titratable groups on mainchain N and H chemical shifts'
        print
        print 'Usage: get_dCS.py <pdbfile> <method>'
        print 'Method is APBS, PBE (Chresten) or Coulomb'
        print
        os._exit(0)
        
