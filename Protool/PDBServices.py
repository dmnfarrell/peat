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

class PDBServices:

    def __init__(self):
        try:
            import SOAPpy
            self.server=SOAPpy.SOAPProxy("http://www.pdb.org/pdb/services/pdbws")
        except:
            print 'SOAP services not available'
        try:
            import Bio
        except:
            print 'BioPython services not available'
        return

    def FASTa(self,pdbID,chainID):
        return self.server.fastaStructureIdQuery(pdbID,chainID,1E-50)

    def getChains(self,PDBID):
        return self.server.getChains(PDBID)

    def getPDB(self,pdbID):
        """Get a PDB file"""
        import urllib
        data = urllib.urlopen('http://www.pdb.org/pdb/files/%s.pdb' %pdbID.lower()).read()
        return data.split('\n')
        