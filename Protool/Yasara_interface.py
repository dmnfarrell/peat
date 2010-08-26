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

import sys, os

class yasara:

    def __init__(self,mode=None):
        yasaradir='/Users/nielsen/Desktop/YASARA.app/yasara' # Replace this with your own yasara directory
        if not os.path.isdir(yasaradir):
            yasaradir=os.path.split(yasaradir)[0]
            dirname=os.path.split(yasaradir)[1]
            if dirname.lower()=='yasara.app':
                yasaradir=os.path.join(yasaradir,'yasara')
        sys.path.append(os.path.join(yasaradir,'pym'))
        sys.path.append(os.path.join(yasaradir,'plg'))
        import yasaramodule as yasara
        self.Y=yasara
        if mode=='txt':
            self.Y.info.mode='txt'
            self.Y.Console('Off')
        return
        
    def shortMD(self,pdblines):
        """Do a short MD in vacou to get an ensemble"""
        objs=self.Y.BuildPDB(pdblines)
        self.Y.DelRes('WAT')
        self.Y.DelRes('HOH')
        self.Y.Clean()
        self.Y.CellAuto(extension=10)
        self.Y.ForceField('Amber99')
        self.Y.Boundary('Periodic')
        self.Y.TempCtrl('SteepDes')
        self.Y.Sim('On')
        self.Y.Wait('SpeedMax<1000')
        self.Y.Sim('Off')
        self.Y.TempCtrl('Rescale')
        self.Y.Sim('On')
        raw_input('yes?')
        return
        #self.Y.Temp(298.15,'Yes')
        
if __name__=='__main__':
    Y=yasara()
    fd=open('2LZT.pdb')
    pdblines=fd.readlines()
    fd.close()
    Y.shortMD(pdblines)
        
        
