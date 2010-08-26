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

class pI:

    def calc_pI(self):
        """Calculate the pI for the calculation currently loaded"""
        calc=self.current_selected_calc
        #
        # Get the pH values
        #
        groups=self.calcs[calc]['titcurv'].keys()
        pHs=self.calcs[calc]['titcurv'][groups[0]].keys()
        pHs.sort()
        #
        # Add all titration curves
        #
        intpkas=[]
        sub_pka={}
        for pH in pHs:
            crg=0.0
            for group in groups:
                crg=crg+self.calcs[calc]['titcurv'][group][pH]
            print '%5.2f: %6.3f' %(pH,crg)
        return
