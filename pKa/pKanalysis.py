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

class pKaIO:

    def __init__(self):
        return

    def readtitcurv(self,filename):
        #
        # This function reads a WHAT IF titration curve file and
        # creates self.titdata, which is a dictionary:
        # self.titdata={<residue>:{'pKa':<pka>,<ph1>:<charge1>,<ph2>:<charge2>.....}}
        #
        import os, string
        if not os.path.isfile(filename):
            raise 'File does not exist:',filename
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        #
        # Parse
        #
        if string.lower(string.strip(lines[0]))!=string.lower('WHAT IF Titration Curve File'):
            raise 'Not a WHAT IF Titration Curve File: ',filename
        if string.lower(string.strip(lines[1]))!=string.lower('Format 1.0'):
            raise 'unknown format: ',lines[1]
        phvals=string.split(lines[2])
        phstart=string.atof(phvals[0])
        phend=string.atof(phvals[1])
        phstep=string.atof(phvals[2])
        titdata={}
        linenumber=3
        done=0
        while not done:
            residue=string.split(lines[linenumber])
            if string.lower(residue[0])==string.lower('TERMINAL'):
                linenumber=linenumber+1
                residue=string.split(lines[linenumber])
                pKa=string.atof(residue[-1])
                residue='T'+string.zfill(residue[0],4)+residue[1]
            else:
                pKa=string.atof(residue[-1])
                residue=string.zfill(residue[0],4)+residue[1]
            linenumber=linenumber+1
            #print 'Reading for: ',residue
            charge={'pKa':pKa}
            for pH in range(100*phstart,100*phend,100*phstep):
                rpH=float(pH)/100.0
                line=string.split(lines[linenumber])
                if string.atof(line[0])==rpH:
                    charge[rpH]=string.atof(line[1])
                    linenumber=linenumber+1
            titdata[residue]=charge
            linenumber=linenumber+2
            if string.strip(string.lower(lines[linenumber]))==string.lower('End of file'):
                done=1
        self.titdata=titdata
        return self.titdata


class pKanalyse:

    def __init__(self):
        return
            
    def calcpKa(self,pH,charge):
        import math
        return pH+math.log(abs(charge)/(1.0-abs(charge)))

    def anacurve(self,residue):
        charges=self.titdata[residue]
        phlist=charges.keys()
        phlist.sort()
        phlist=phlist[:-1]
        pKa={}
        list=[]
        #print
        #print 'Residue: ',residue
        #print 'pKa given by WHAT IF: %6.3f' %(charges['pKa'])
        #print
        #print '   pH    pKa   charge '
        for pH in phlist:
            charge=charges[pH]
            if abs(charge)>0.20 and abs(charge)<0.80:
                pKa[pH]=self.calcpKa(pH,charge)
                list.append(pKa[pH])
                #print ' %4.2f  %6.3f %5.3f' %(pH, pKa[pH],charges[pH])
            else:
                pKa[pH]='****'
        #print
        avg,max,min=self.stats(list)
        return max-min

    def stats(self,list):
        sum=0.0
        max=-9999
        min=9999
        if len(list)!=0:
            for value in list:
                sum=sum+value
                if value>max:
                    max=value
                if value<min:
                    min=value
            return sum/len(list),max,min
        else:
            return 0.0,max,min

    

    
    
