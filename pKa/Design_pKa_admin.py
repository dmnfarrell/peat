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

def get_mut_name(mutation,Protool_instance):
    """Get a mutation name - helper function for Design_pKa reading of WHAT IF SUGELM file
    Protool_instance must have read the original PDB file
    Produces a mutation name of teh type A:0010:ASP:ASN for the mutatino of ASP 10 in chain A to ASN"""
    #
    # Get the residue number (includes chain ID)
    #
    import pKD_tools
    resnum=pKD_tools.get_resid_from_mut(mutation)
    #
    # Get the new residue type
    #
    newres=mutation.split(':')[-1]
    #
    # New name
    #
    newname='%s:%s:%s' %(resnum,Protool_instance.resname(resnum),newres)
    return newname
