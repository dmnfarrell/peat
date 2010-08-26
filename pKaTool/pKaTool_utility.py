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


class pKaTool_utility:

    def get_charge(self,group,target_pH,calc=None):
        """Get the charge of a specific group at a specific pH.
        If the data point is not there, then do a linear interpolation"""
        if not calc:
            calc=self.current_selected_calc
        if not calc:
            print 'No calc selected'
            return None
        #
        # Find the closest pHs for the group
        #
        charges=self.calcs[calc]['titcurv'].copy()
        best_pHs=[]
        for pH in charges[group].keys():
            if pH!='pKa':
                best_pHs.append([abs(float(pH)-target_pH),pH,charges[group][pH]])
        best_pHs.sort()
        #
        # Do we have an exact match?
        #
        if abs(best_pHs[0][1]-target_pH)<0.05:
            return float(best_pHs[0][2])
        else:
            #
            # No, do interpolation between two best values
            # but only if they both are closer than 1.0 unit
            #
            if abs(best_pHs[0][1]-target_pH)<1.0 and abs(best_pHs[1][1]-target_pH)<1.0:
                charge_diff=best_pHs[1][2]-best_pHs[0][2]
                pH_diff=best_pHs[1][1]-best_pHs[0][1]
                slope=charge_diff/pH_diff
                guessed_charge=(target_pH-best_pHs[0][1])*slope+best_pHs[0][2]
                return guessed_charge
            else:
                print 'Target: %5.2f  pH1: %5.2f  pH2: %5.2f - returning none' %(target_pH,best_pHs[0][1],best_pHs[1][1])
                return None

    

                             
