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


class titration_curve:

    def __init__(self,curves):
        self.curves=curves.copy()
        not_allowed=['pKa','pka']
        for key in not_allowed:
            if self.curves.has_key(key):
                del self.curves[key]
        return

    #
    # ----
    #

    def __sub__(self,other):
        """subtract two titration curves"""
        diff=0.0
   
        for group in self.curves.keys():
            if not other.curves.has_key(group):
                continue
            for ph in self.curves[group].keys():
                if other.curves[group].has_key(ph):
                    diff=diff+abs(self.curves[group][ph]-other.curves[group][ph])
        return diff

    #
    # ----
    #

    def subtract_individually(self,other):
        """Subtract curves individually"""
        diff=[] 
        for group in self.curves.keys():
            if not other.curves.has_key(group):
                diff.append(0.0)
                continue
            this_diff = 0.0
            for ph in self.curves[group].keys():
                if other.curves[group].has_key(ph):
                    this_diff=this_diff+abs(self.curves[group][ph]-other.curves[group][ph])
            diff.append(this_diff)
        return diff

    #
    # ----
    #

    def sub_scaled(self,other):
        """scaled difference btw two titration curves"""
        diff=0.0
        for group in self.curves.keys():
            if not other.curves.has_key(group):
                continue
                raise 'incompatible titration curves'
            for ph in self.curves[group].keys():
                if other.curves[group].has_key(ph):
                    diff=diff+self.scale(self.curves[group][ph],other.curves[group][ph])*abs(self.curves[group][ph]-other.curves[group][ph])
        return diff

    #
    # ----
    #

    def scale(self,frac1,frac2):
        """Scales the error on a titration point """
        return max(self.scale_function(frac1),self.scale_function(frac2))

    #
    # ----
    #

    def scale_function(self, x):
        """Calculates the scaling functionfor scaling """
        return -pow(abs(x)-1,2)+1

    #
    # ----
    #

    def experimental_uncertainty(self, pH_uncertainty=0.1):
        """estimates the experimental uncertainty of titration curves"""
        res=0.0
        count = 0
        for group in self.curves.keys():
            #print 'Now estimating for ',group
            pHs = self.curves[group].keys()
            #make sure that ph values are sorted
            pHs.sort()
            for i in range(len(pHs)):
                bw_diff = 0
                fw_diff = 0
                try:
                    bw_diff = (self.curves[group][pHs[i]]-self.curves[group][pHs[i-1]])/(pHs[i]-pHs[i-1])
                except:
                    pass
                try:
                    fw_diff = (self.curves[group][pHs[i+1]]-self.curves[group][pHs[i]])/(pHs[i+1]-pHs[i])
                except:
                    pass

                avr_diff = (bw_diff+fw_diff)/2 ##### abs()?
                res += avr_diff
                count += 1
        res *= pH_uncertainty
        res = abs(res)
        avr_res = res / float(count)
        return res, avr_res

    #
    # ----
    #

    def sub_HHd_scaled(self, exp_data, pkas):
        """Calculates error with scaling based on deviation of exp data from the Henderson-Hasselbalch eq"""
        diff=0.0
        scales = exp_data.deviation_from_henderson_hasselbalch(pkas)
        for group in self.curves.keys():
            if not exp_data.curves.has_key(group):
                continue
            for ph in self.curves[group].keys():
                if exp_data.curves[group].has_key(ph):
                    diff=diff+ scales[group][ph]*abs(self.curves[group][ph]-exp_data.curves[group][ph])
        return diff

    #
    # -----
    #

    def deviation_from_henderson_hasselbalch(self, pKas):
        """Calculates the deviation from the Henderson-Hasselbalch equation for all points given pKa values"""
        HH_deviation = {}
        deviation = lambda ph,pka,exp: abs(1/(1+pow(10,ph-pka))-1-exp)
        for group in self.curves.keys():
            if pKas.has_key(group):
                pka = pKas[group]
                HH_deviation[group] = {}
                for ph in self.curves[group].keys():
                    try:
                        HH_deviation[group][ph] = deviation(float(ph),float(pka),float(self.curves[group][ph]))
                    except:
                        pass
        return HH_deviation
    
