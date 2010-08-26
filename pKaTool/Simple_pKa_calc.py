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

import pKaTool.pKa_calc as pKa_calc

def main():
    X=pKa_calc.Boltzmann() # You can use other classes instead of Boltzmann if you want to use a different method for determinig the titration curves:
    # Tanford_Roxby
    # Monte_Carlo
    # Monte_Carlo_CPP (C++ implementation of Monte Carlo routine - much faster)
    # Boltzmann_CPP (C++ implementation of the Boltzmann routine)
    
    #
    # Fill in desolvation energies and site-site interaction energies
    #
    X.groups=[':0001:ACID',':0002:BASE']
    X.intene={}
    X.intrinsic_pKa={}
    for group in X.groups:
        import random
        X.intrinsic_pKa[group]=random.uniform(4,10) # Insert intrinsic pKa values here
        X.intene[group]={}
        for group2 in X.groups:
            X.intene[group][group2]=0.0 # Insert interaction energies here - remember that matrix must be symmetric
    pKa_values,protstates=X._calc_pKas()
    print 'pKa values',pKa_values
    #
    # Plot the titration curves using matplotlib
    #
    try:
        import pylab
    except:
        return
    for group in X.groups:
        pHs=protstates[group].keys()
        pHs.sort()
        pHs=pHs[:-1]
        charge=[]
        for pH in pHs:
            charge.append(protstates[group][pH])
        pylab.plot(pHs,charge,label=group)
    pylab.legend()
    pylab.show()
    return
    

if __name__=='__main__':
    main()