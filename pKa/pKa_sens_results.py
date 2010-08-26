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
import string, math, os, access, math, pKa, ostools

def avg_charge_pKa(files):
    #
    # Calculate the average titration curve and the average pKa for the dirs
    #
    import sys, os, string
    #dirs=sys.argv[1:]
    top=os.getcwd()
    pkadict={}
    titdict={}
    ok=0
    wrong=0
    for file in files:
        if os.path.isfile(file):
            #print 'Processing file:',file
            import pKa
            if os.path.isfile(file):
                X=pKa.pKaIO()
                tit=X.readtitcurv(file)
                titdict[file]={}
                if tit.has_key('0035GLU') and tit.has_key('0052ASP'):
                    titdict[file]={'35':tit['0035GLU'],'52':tit['0052ASP']}
                    
            else:
                print 'pKa calculations not finished for: %s' %dir
    #
    # Calculate averages
    #
    sumtit={}
    #
    # Zero dictionaries
    #
    for outkey in titdict.keys():
        for key in titdict[outkey].keys():
            sumtit[key]={}
    #
    # Sum all
    #
    pkacount=0
    for dir in titdict.keys():
        for key in titdict[dir].keys():
            for ph in titdict[dir][key].keys():
                if sumtit[key].has_key(ph):
                    sumtit[key][ph]=sumtit[key][ph]+titdict[dir][key][ph]
                else:
                    sumtit[key][ph]=titdict[dir][key][ph]
    #
    # Get average 
    #
    for key in sumtit.keys():
        for ph in sumtit[key].keys():
            sumtit[key][ph]=float(sumtit[key][ph])/float(len(titdict.keys()))
    #
    # Print the results
    #
    residues=sumtit.keys()
    residues.sort()
    pka=0.0
    print 'pKa value from average titration curve'
    for key in residues:
        phvals=sumtit[key].keys()
        phvals.sort()
        for ph in phvals:
            if sumtit[key][ph]<=-0.5:
                pka=ph
                break
        print "Residue: %s, pKa value: %5.2f" %(key,pka)
    return

def do_dir(path,detailed=None):
    #
    # Load experimental pKa values
    #
    expdata='/u1/jnielsen/work/pka/exppka/lysozyme.PKA.DAT'
    experimental=pKa.pKaIO()
    experimental.readpka(expdata)
    A=pKa.pKanalyse()
    #
    # Find all directories with pKa calculations and loop over them
    #
    list=ostools.find('*.PKA.DAT',path)
    topdir=os.getcwd()
    print
    print '================================================'
    print 'Evaluating ',len(list),' calculations in: %s' %os.path.split(path)[1]
    #
    # Reset counters
    #
    correct=0
    wrong=0
    rmsds=[]
    E35sum=[]
    D52sum=[]
    for file in list:
        dir,pkafile=os.path.split(file)
        pos=string.find(dir,path)
        os.chdir(dir)
        if pos!=1:
            dir=dir[pos+len(path):]
        calculated=pKa.pKaIO()
        calculated.readpka(pkafile)
        rmsd=A.rmsd_between_sets(experimental.pka,calculated.pka)
        if detailed:
            print 'Dir: %6s rmsd: %5.2f ' %(dir,rmsd)
            #
            # Get stats on D52 and E35
            #
            residues=['0035GLU','0052ASP']
            print '          Residue  Exp. Calc. diff'
            for residue in residues:
                print ' %15s %5.2f %5.2f %5.2f' %(residue,experimental.pka[residue]['pKa'], calculated.pka[residue]['pKa'],experimental.pka[residue]['pKa']-calculated.pka[residue]['pKa'])
            if calculated.pka['0035GLU']['pKa']-calculated.pka['0052ASP']['pKa']>=1.5 and calculated.pka['0035GLU']['pKa']>5.0:
                print 'GLU35!'
            elif -calculated.pka['0035GLU']['pKa']+calculated.pka['0052ASP']['pKa']>=1.5 and calculated.pka['0052ASP']['pKa']>5.0:
                print 'ASP52!'
            else:
                print 'None'
            print
        #
        # Did we get the h-donor right?
        #
        if calculated.pka.has_key('0035GLU'):
            E35=calculated.pka['0035GLU']['pKa']
            E35sum.append(E35)
        if calculated.pka.has_key('0052ASP'):
            D52=calculated.pka['0052ASP']['pKa']
            D52sum.append(D52)
        if calculated.pka.has_key('0035GLU') and calculated.pka.has_key('0052ASP'):
            if E35>5.0 and E35-D52>=1.5:
                correct=correct+1
            if D52>5.0 and D52-E35>=1.5:
                wrong=wrong+1
        rmsds.append(rmsd)
    #
    # --------------------------------------------------------
    #
    # Print the summary
    #
    print 'Number of dirs',len(rmsds)
    if len(rmsds)>0:
        sum=0.0
        for rmsd in rmsds:
            sum=sum+rmsd
        Esum=0.0
        for p in E35sum:
            Esum=Esum+p
        Dsum=0.0
        for p in D52sum:
            Dsum=Dsum+p
        D52mean=Dsum/float(len(D52sum))
        E35mean=Esum/float(len(E35sum))
        avgdevE=0.0
        avgdevD=0.0
        for p in E35sum:
            avgdevE=avgdevE+abs(p-E35mean)
        for p in D52sum:
            avgdevD=avgdevD+abs(p-D52mean)
        avgdevE=avgdevE/float(len(E35sum))
        avgdevD=avgdevD/float(len(D52sum))
        print 'Average rmsd: %5.2f' %(sum/float(len(rmsds)))
        print 'Average pKa E35: %5.2f' %(Esum/float(len(E35sum)))
        print 'Average pKa D52: %5.2f' %(Dsum/float(len(D52sum)))
        print 'Avg. deviation from mean (E35): %7.3f' %avgdevE
        print 'Avg. deviation from mean (D52): %7.3f' %avgdevD
        print 'Correct:',correct
        print 'Wrong:',wrong

        #
        # Average the titration curves
        #
        list=ostools.find('*.TITCURV.DAT',path)
        avg_charge_pKa(list)
        os.chdir(topdir)
    return



def main():
    #
    # Run for several directories
    #
    import sys
    detailed=None
    if len(sys.argv)>1:
        detailed=1
    dirs=['Xray','EM','MD','MD_structs_EM','averagestructs','averagestructs_EM','long_MD/2LZT','long_MD/2LZT_D52H','long_MD/2LZT_E35H','long_MD/2LZT_E35H_D52H','long_MD/5LYZ','long_MD/2LZT_reallong','Xray_indi4','Xray_indi16','NMR','random']
    
    dirs=['Xray_indi16']
    #dirs=['Xray']
    cwd=os.getcwd()
    for dir in dirs:
        dir=os.path.join(cwd,dir)
        do_dir(dir,detailed)
    return

main()
