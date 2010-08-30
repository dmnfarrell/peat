#!/usr/bin/env python
#
# FFF - Flexible Force Field
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


import pylab, math
k=1.3806503E-23
T=298.15
Na=6.02E23
kT=k*T*Na/1000.0
print '1 kT is %6.2f kJ/mol' %kT
print '1 dpKa is %6.2f kJ/mol' %(kT*math.log(10))
class Plots:

    def __init__(self,options):
        """Make the plots specified"""
        self.pkacalcerror=0.7
        self.IEerr=0.2
        #
        options=eval((str(options)))
        for func,numpoints,value in [['pkavalues',108,'N']]:#,['pkavalues',250,'H'],['pkavalues',100,'C'],['intenes',10,None],['dpkas',30,None]]:
            if not value:
                xs,ys,col,legend=getattr(self,func)(numpoints)
            else:
                xs,ys,col,legend=getattr(self,func)(numpoints,value)
                pylab.plot(xs,ys,col,label=legend)
        pylab.xlabel('Experimental values')
        pylab.ylabel('Calculated values')
        pylab.legend()
        pylab.show()
        return
        
    def pkavalues(self,numpoints,nucleus='N'):
        """Make a plot of experimental values including their estimated errors"""
        #
        # pKa value shifts range from <-4 -> +4 with an error of 0.5 units
        #
        import numpy
        xs=[]
        ys=[]
        error={'N':0.5,'C':0.1,'HN':0.5,'H':0.2}
        cols={'N':'b','C':'c','H':'g','HN':'k'}
        for count in range(numpoints):
            import random
            dpka=random.normalvariate(0.0,0.5)
            dpkay=dpka+random.normalvariate(0.0, error[nucleus])
            calpka=dpka+random.normalvariate(0.0,self.pkacalcerror)
            ys.append(dpka*math.log(10)*kT)
            xs.append(dpkay*math.log(10)*kT)
        return xs,ys,'%1so' %cols[nucleus],'pKas (%s)' %nucleus
        
    def intenes(self,numpoints):
        """Electrostatic interaction energies from GloFTE fits"""
        import numpy
        xs=[]
        ys=[]
        for count in range(numpoints):
            import random
            intene=random.normalvariate(0.0,1.0)
            inteney=intene+random.normalvariate(0.0, 0.25)
            ys.append(intene*kT)
            xs.append(inteney*kT)
        return xs,ys,'bo','Eints'

    def dpkas(self,numpoints):
        """Electrostatic interaction energies from GloFTE fits"""
        import numpy
        xs=[]
        ys=[]
        for count in range(numpoints):
            import random
            intene=random.normalvariate(0.0,0.3)
            inteney=intene+random.normalvariate(0.0, 0.05)
            ys.append(intene*math.log(10)*kT)
            xs.append(inteney*math.log(10)*kT)
        return xs,ys,'go','dpKas'


if __name__=='__main__':
    from optparse import OptionParser
    import os
    parser = OptionParser(usage='%prog [options]',version='%prog 1.0')
    (options, args) = parser.parse_args()
    X=Plots(options)
