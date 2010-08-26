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

from sys import argv
import string, math

def getdiff(pkafile,exppka=None):
    if not pkafile:
        raise "Too few arguments"
    if exppka:
        file1=exppka
        file2=pkafile
    else:
        if len(argv)>2:
            file1=pkafile
            file2=argv[2]
        else:
            file1='/u1/jnielsen/work/pka/exppka/lysozyme.PKA.DAT'
            file2=pkafile
    """Read the two files and compute RMS"""
    fd1=open(file1)
    fd2=open(file2)
    lin1=fd1.readlines()
    lin2=fd2.readlines()
    if len(lin1)!=len(lin2):
        print len(lin1),len(lin2)
        print 'Files of unequal length'
        return -1, -1, -1
    else:
        print
        print 'Exp. values',file2,'Difference'
        rmsd=0.0
        rmsd2=0.0
        count=0.0
        missing=0
        htmlstring=''
        diffstring=''
        for line in range(len(lin1)-2):
            word1=string.split(lin1[line+2])
            word2=string.split(lin2[line+2])
            if word1[:5]==word2[:5] and string.strip(word1[0])[0] in string.digits:
                val1=word1[5][0]
                val2=word2[5][0]
                #print val1, ':',val2
                if val1 not in string.digits:
                    val1=word1[6]
                if val2 not in string.digits:
                    val2=word2[6]
                diff=string.atof(val1)-string.atof(val2)
                pka1=string.atof(val1)
                pka2=string.atof(val2)

                    
                if pka2==0.0:
                    diff=0
                    if string.strip(word1[1])!='ARG':                
                        missing=missing+1
                if (pka1==0):
                    diff=0.0
                if string.strip(word1[0])=='35' or string.strip(word1[0])=='52':
                    rmsd2=rmsd2+math.pow(diff,2)
                if pka2!=0.0:
                    text="%15s %5.2f %5.2f %6.2f" %(string.join(word1[:5]),pka1,pka2,diff)
                    htmlstring=htmlstring+"<TD> %5.2f" %pka2
                    diffstring=diffstring+"<TD> %6.2f" %diff
                else:
                    text="%15s %5.2f %5s %6s" %(string.join(word1[:5]),pka1,'****','****')
                    htmlstring=htmlstring+"<TD> %6s" %('****')
                    diffstring=diffstring+"<TD> %6s" %('****')
                print text
                if (pka1!=0.0 and pka2!=0.0):
                    rmsd=rmsd+math.pow(diff,2)
                    count=count+1.0
        print
        print 'rmsd: %5.2f' %(math.sqrt(rmsd/count))
        print 'rmsd for E35 and D52 %5.2f' %(math.sqrt(rmsd2/count))
        print 'pKa values outside range: ',missing
        print
        htmlstring=htmlstring+'<TD BGCOLOR="pink">%5.2f' %(math.sqrt(rmsd/count))
        htmlstring=htmlstring+'<TD>%5.2f' %(math.sqrt(rmsd2/count))
        htmlstring=htmlstring+'<TD>'+str(missing)
        diffstring=diffstring+'<TD BGCOLOR="pink">%5.2f' %(math.sqrt(rmsd/count))
        diffstring=diffstring+'<TD>%5.2f' %(math.sqrt(rmsd2/count))
        diffstring=diffstring+'<TD>'+str(missing)
        print htmlstring
        print
        print 'Differences:'
        print diffstring
        return (math.sqrt(rmsd/count)),(math.sqrt(rmsd2/count)),missing
            

if __name__ == "__main__":
    rmsd1,rmsd2,miss=getdiff(argv[1])
