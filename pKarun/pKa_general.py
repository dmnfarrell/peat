#!/usr/bin/env python
#
# pKarun - scripts for running pKa calculations
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

from pKarun.pKa_utility_functions import *
Error='pKa Error'

from numpy import *
from numpy.linalg import *
linear_least_squares=lstsq


class pKanalyse:

    def __init__(self):
        import pKaTool.pKaIO
        X=pKaTool.pKaIO.pKaIO()
        self.modelpk=X.modelpKas
        self.acidbase=X.acidbase
        return

    #
    # -----------
    #
    
    def calcpKa(self,pH,charge):
        import math
        return pH+math.log(abs(charge)/(1.0-abs(charge)))

    #
    # -----------
    #

    def calcpKa_from_titration_curve(self,dict,halfcharge=0.5):
        #
        # Find the pKa value for a titration curve given in dict
        # dict has the format: {ph1:charge1,ph2:charge2,...}
        #
        import string
        phvals=dict.keys()
        phvals.sort()
        if string.lower(str(phvals[-1]))=='pka':
            phvals=phvals[:-1]
        lastcha=abs(dict[phvals[0]])
        pKa=None
        phstep=phvals[1]-phvals[0]
        for ph in phvals[1:]:
            nowch=abs(dict[ph])
            if nowch>=halfcharge and lastcha<=halfcharge and (nowch-lastcha)>0.0:
                pKa=((halfcharge-lastcha)/(nowch-lastcha))*phstep+(ph-phstep)
            if nowch<=halfcharge and lastcha>=halfcharge and (lastcha-nowch)>0.0:
                pKa=((halfcharge-nowch)/(lastcha-nowch))*phstep+(ph-phstep)
            lastcha=nowch
        #
        # -----------------
        #
        if not pKa:
            pKa=99.0
        return pKa
    
    #
    # ----------------
    #

    def calc_dpka(self,dict,subset=None):
        #
        # Calculates the delta pKa value from the pKa values in dict
        #
        # Subset can be used to specify a subset of residues to work with
        #
        dpkas={}
        residues=dict.keys()
        if subset:
            residues=subset
        for residue in residues:
            import types
            if type(residue) is types.ListType:
                #
                # If it is a list of residues, then we must first get the pKa value
                #
                X=pKagraph(dict)
                adddict,status=X.addtitcurvs(residue)
                if status:
                    #
                    # Work out what the half-point of the titration should be
                    #
                    halfpoint=0.0
                    for res in residue:
                        if res[-3:]=='ASP' or res[-3:]=='GLU' or res[-3:]=='TYR' or res[-3:]=='CYS':
                            halfpoint=halfpoint-0.5
                        else:
                            halfpoint=halfpoint+0.5
                    halfpoint=abs(halfpoint)
                    #
                    # Special for the chitinase calcs:
                    #
                    halfpoint=-1.5
                    pka=self.calcpKa_from_titration_curve(adddict[str(residue)],halfpoint)
                    dpkas[str(residue)]={'pKa':pka,'dpKa':None}
                else:
                    dpkas[str(residue)]={'pKa':None,'dpKa':None}
            elif not dict.has_key(residue):
                dpkas[residue]={'pKa':None,'dpKa':None}
                continue
            elif residue[0]!='T':
                mpka=self.modelpk[residue[-3:]]
                dpkas[residue]={'pKa':dict[residue]['pKa'],'dpKa':dict[residue]['pKa']-mpka[0]}
            else:
                #
                # Work out if it is an acid or a base
                #
                acid=None
                for pH in dict[residue].keys():
                    if dict[residue][pH]<0.0:
                        acid=1
                #
                # If it is an acid then we guess that it is a C-terminal
                #
                term='NT'
                if acid:
                    term='CT'
                dpkas[residue]={'pKa':dict[residue]['pKa'],'dpKa':dict[residue]['pKa']-self.modelpk[term][0]}
                
        return dpkas

    #
    # ------------
    #
    
    def anacurve(self,residue):
        charges=self.titdata[residue]
        phlist=charges.keys()
        phlist.sort()
        phlist=phlist[:-1]
        pKa={}
        list=[]
        for pH in phlist:
            charge=charges[pH]
            if abs(charge)>0.20 and abs(charge)<0.80:
                pKa[pH]=self.calcpKa(pH,charge)
                list.append(pKa[pH])
            else:
                pKa[pH]='****'
        avg,max,min=self.stats(list)
        return max-min

    #
    # -----------------
    #

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

    #
    # ------
    #

    def fit_to_henderson(self,curve):
        """Fit a titration curve to the Henderson-Hasselbalch equation: pH = b + a(log([HA]/[A-]))
        a (the slope) will be 1 for a perfect HH curve, and b will be the pKa value"""
        import string, math, os

        keys=curve.keys()
        keys.sort()
        if string.strip(string.lower(str(keys[-1])))=='pka':
            keys=keys[:-1]
        #
        # curve is a dictionary of pH values:charge
        #
        matrix=[]
        vector=[]
        weights={}
        logcurve={}
        for pH in keys:
            if abs(curve[pH])>0.05 and abs(curve[pH])<0.95:
                logcurve[pH]=math.log10(1.0/abs(curve[pH])-1.0)
                weights[pH]=1.0
        for pH in logcurve.keys():
            matrix.append([pH,weights[pH]])
            vector.append(logcurve[pH])
        matrix=array((matrix))
        vector=array((vector))
        if len(logcurve.keys())<2:
            return [0.0,0.0],[0.0]
        try:
            solution,sq,rank,xx=linear_least_squares(matrix,vector)
        except "LinearAlgebraError":
            print 'ERROR', logcurve.keys()
            print matrix,vector
            solution=[0.0,0.0]
            sq=[0.0]
        return solution,sq

    #
    # ----------
    #

    def fit_all_to_henderson(self,dict,subset=None):
        #
        # dict is a dictionary that ???
        # Fits all titration curves in dict to a Henderson-Hasselbalch curve
        #
        data={}
        keys=dict.keys()
        keys.sort()
        if subset:
            keys=subset
        for residue in keys:
            import types
            if type(residue) is types.ListType:
                X=pKagraph(dict)
                adddict,status=X.addtitcurvs(residue)
                if status:
                    solution,sq=self.fit_to_henderson(adddict[str(residue)])
                    try:
                        data[str(residue)]=[solution[0],solution[1],sq[0]]
                    except IndexError or TypeError:
                        data[str(residue)]=[solution[0],solution[1],0]
                else:
                    data[str(residue)]=[None,None,None]
            elif dict.has_key(residue):
                solution,sq=self.fit_to_henderson(dict[residue])
                try:
                    data[residue]=[solution[0],solution[1],sq[0]]
                except IndexError or TypeError:
                    data[residue]=[solution[0],solution[1],0]
            else:
                data[residue]=[None,None,None]
        return data

    #
    # -----
    #
            
    def analyse_stability(self,folded,unfolded=0,absolute=False,write_file=True):
        """
        # folded is a dictionary containing the titration
        # curves for all residues in the protein. (Normally
        # the self.titdata produced by pKaIO.readtitcurv).
        # unfolded is a similar dictionary, containing the
        # titration curves for the unfolded form. If no unfolded
        # dictionary is given, then titration curves
        # are calculated based on the model pKa values for
        # the amino acid types.
         """
        import string, math
        residues=folded.keys()
        checkedres=[]
        for res in residues:
            if string.strip(string.lower(res))!='pka':
                checkedres.append(res)
        residues=checkedres
        if unfolded==0 or not unfolded:
            unfolded={}
            X=pKatools()
            for residue in residues:
                type=self.get_residue_type(residue)
                
                phstart,phend,phstep,pka,curve=self.get_residue_curve(folded[residue])
                unfolded[residue]=X.make_titration_curve_2(self.modelpk[string.upper(type)],self.acidbase[string.upper(type)],residue,folded[residue].keys())[residue]
        new_phvalues=folded[residues[0]].keys()
        phvalues=[]
        for ph in new_phvalues:
            if string.lower(str(ph))!='pka' and ph>=2.0 and ph<=10.0:
                phvalues.append(ph)
        phvalues.sort()
        stab={}
        RK=8.31451
        TEMP=298.15
        #
        # Get the stability "titration" curve for each residue
        #
        sums={}
        for residue in residues:
            sumc=0.0
            lastph=phvalues[0]
            contribution={}
            for ph in phvalues[1:]:
                print ph,folded[residue][ph],unfolded[residue][ph],residue
                diff=math.log(10)*RK*TEMP*(folded[residue][ph]-unfolded[residue][ph])/1000.0 #Result in kJ/mol?
                if absolute:
                    diff=abs(diff)
                if self.is_base(residue):
                    diff=-diff
                sumc=sumc+diff*(ph-lastph)
                contribution[ph]=sumc
                lastph=ph
            sums[residue]=sumc
            stab[residue]=contribution
        # Done
        return stab,sums

    def get_residue_curve(self,curve):
        import string
        keys=curve.keys()
        keys.sort()
        pka=-1
        for key in keys:
            if string.strip(string.lower(str(key)))=='pka':
                pka=curve[key]
                del curve[key]
        phvals=curve.keys()
        phvals.sort()
        phstart=phvals[0]
        phend=phvals[-1]
        phstep=abs(phvals[1]-phvals[0])
        return phstart,phend,phstep,pka,curve

    #
    # ----
    #

    def get_residue_type(self,residue):
        """From the residue name, get the residue type"""
        rtype=residue.split(':')[-1]
        return rtype

    #
    # ----
    #

    def is_acid(self,residue):
        type=self.get_residue_type(residue)
        if self.acidbase[type]==-1:
            return True
        return False

    def is_base(self,residue):
        type=self.get_residue_type(residue)
        if self.acidbase[type]==1:
            return True
        return False    

    def rmsd(self,list):
        import math
        rmsd=0.0
        for value in list:
            rmsd=rmsd+math.pow(value,2.0)
        return math.sqrt(rmsd/len(list))

    def rmsd_between_sets(self,calc1,calc2):
        #
        # This function calculates the rmsd between two sets
        # of pKa values. calc1 and calc2 are dictionaries as
        # defined by readpka in pKaIO
        # Only residues common to both sets are compared!
        # Residues found in one set, but not the other, are reported
        #
        # If any of the pKa values are zero, then the difference for this
        # group is not reported
        #
        diflist=[]
        zero1=0
        zero2=0
        keys1=calc1.keys()
        keys1.sort()
        for residue in keys1:
            if calc2.has_key(residue):
                if calc1[residue]['pKa']!=0.0 and calc2[residue]['pKa']!=0.0:
                    diflist.append(calc1[residue]['pKa']-calc2[residue]['pKa'])
            else:
                pass
            if int(calc1[residue]['pKa'])==0:
                zero1=zero1+1
        for residue in calc2.keys():
            if int(calc2[residue]['pKa'])==0:
                zero2=zero2+1
            if not calc1.has_key(residue):
                pass
        self.diflist=diflist
        return self.rmsd(diflist)

    def calculate_pI(self,sequence,nocys=None):
        #
        # Calculates the pI of a protein sequence. The model pKas are taken from the
        # self.modelpK array
        # sequence is a string containing the protein sequence in one-letter code
        #
        phvals=range(000,1500,10)
        lastcharge=+100.0
        for phval in phvals:
            ph=float(phval)/100.0
            charge=0.0
            for residue in sequence:
                if residue=='K':
                    charge=charge+self.getcharge('LYS',ph)
                elif residue=='R':
                    charge=charge+self.getcharge('ARG',ph)
                elif residue=='H':
                    charge=charge+self.getcharge('HIS',ph)
                elif residue=='D':
                    charge=charge+self.getcharge('ASP',ph)
                elif residue=='E':
                    charge=charge+self.getcharge('GLU',ph)
                elif residue=='C':
                    chargec=self.getcharge('CYS',ph)
                    if nocys:
                        chargec=0.0
                    charge=charge+chargec
                elif residue=='Y':
                    charge=charge+self.getcharge('TYR',ph)
                else:
                    pass
            charge=charge+self.getcharge('NT',ph)+self.getcharge('CT',ph)
            if charge<0.0 and lastcharge>=0.0:
                frac=(0.0-lastcharge)/(charge-lastcharge)
                pI=ph-0.1+0.1*frac
                pI='%5.3f' %pI
                return pI
            else:
                lastcharge=charge
        return None

    def getcharge(self,residue,ph):
        #
        # Calculates the charge on a particular group at a
        # certain pH using the Henderson-Hasselbalch equation
        #
        if self.modelpk.has_key(residue):
            pklist=self.modelpk[residue]
            pka=pklist[0]
            if pklist[1]==0:
                #
                # Base
                #
                charge=1.0/(10.0**(ph-pka)+1.0)
                return charge
            else:
                #
                # Acid
                #
                charge=-(1.0-(1.0/(10.0**(ph-pka)+1.0)))
                return charge
        else:
            return None

        
#======================================================================
class pKagraph:

    def __init__(self,tdata):
        self.tdata=tdata
        return

    def addtitcurvs(self,reslist):
        #
        # Add two titration curves and return a new dictionary
        #
        # reslist is a list e.g. ['0138ASP','0140ASP']
        #
        print 'Adding titration curves for',reslist
        adddict={str(reslist):{}}
        for res in reslist:
            if not self.tdata.has_key(res):
                return '@'+str(reslist),None
            for ph in self.tdata[res].keys():
                if not adddict[str(reslist)].has_key(ph):
                    adddict[str(reslist)][ph]=0.0
                adddict[str(reslist)][ph]=adddict[str(reslist)][ph]+self.tdata[res][ph]
        return adddict,1

    def plotcurv(self,tdata,plot_henderson=1,dir='gifs',subset=None):
        #
        # Tdata is a dictionary {residue:{pH value:charge}}.
        # The subroutine will loop over all residues and create
        # a gif for each of them.
        # The entry 'pKa' is removed automatically from tdata
        # If plot_henderson==1 then will a perfect Henderson-
        # Hasselbalch curve be displayed too.
        # dir controls the subdir that will be created to store the
        # gif files.
        # The routine returns a list of gif filenames.
        #
        print 'Plotting titration curves for residues:'
        import os, sys
        giffiles=[]
        #
        # Loop over all residues
        #
        import types
        residues=tdata.keys()
        if subset:
            residues=subset
        residues.sort()
        for residue in residues:
            #
            # If we are asked to add some titration curves
            #
            if type(residue) is types.ListType:
                adddict,status=self.addtitcurvs(residue)
                if status:
                    giffiles.append(self.plot_single_curve(adddict,str(residue),plot_henderson,dir))
                else:
                    giffiles.append(adddict)
            elif tdata.has_key(residue):
                giffile=self.plot_single_curve(tdata,residue,plot_henderson,dir)
                giffiles.append(giffile)
            else:
                giffiles.append('@'+residue)
        return giffiles


    
    def plot_single_curve(self,tdata,residue=None,plot_henderson=1,dir='gifs'):
        #
        # Plots a single curve in tdata
        # Plots a perfect Henderson-Hasselback curve for comparison
        #
        import sys, Gnuplot, os, string, time
        g=Gnuplot.Gnuplot()
        if residue:
            print residue,
            sys.stdout.flush()
            curve=tdata[residue]
        else:
            curve=tdata
        X=pKanalyse()
        phstart,phend,phstep,pka,curve=X.get_residue_curve(curve)
        plotdata=[]
        max=-99.9
        min=99.9
        phvals=curve.keys()
        phvals.sort()
        for ph in phvals:
            plotdata.append([ph,curve[ph]])
            if curve[ph]>max:
                max=curve[ph]
            if curve[ph]<min:
                min=curve[ph]        
        #data=Gnuplot.Data(plotdata,title=residue,_with='lines')
        #============================================================
        # Set min and max 
        #
        rmax=round(max)
        rmin=round(min)
        if rmax==rmin:
            if max>rmax:
                rmax=rmax+1
            if min<rmin:
                rmin=rmin-1
        if rmax==rmin:
            if residue[-3:]=='ASP' or residue[-3:]=='GLU' or residue[-3:]=='TYR':
                rmin=rmin-1
                rmax=rmax+0.1
            else:
                rmax=rmax+1
                rmin=rmin-0.1
        
        #============================================================
        # Make an ideal Henderson-Hasselbalch curve for comparison
        #
        if pka!=-1 and plot_henderson==1:
            if rmin<0.0 and rmax<=0.0:
                acid=1
            else:
                acid=0
            tools=pKatools()
            dict=tools.make_titration_curve(pka,acid,'ideal',phstart,phend,0.1)
            plotdata2=[]
            keys=dict[dict.keys()[0]].keys()
            keys.sort()
            for ph in keys:
                if string.lower(string.strip(str(ph)))!='pka':
                    plotdata2.append([ph,dict[dict.keys()[0]][ph]])
            #data2=Gnuplot.Data(plotdata2,title='Ideal Henderson-Hasselbalch',_with='lines')
        #
        # Make the plot
        #
        g.title('Titration curve')
        g.xlabel('pH')
        g.ylabel('Charge')
        g('set yrange ['+str(rmin)+':'+str(rmax)+']')
        if pka!=-1 and plot_henderson==1:
            g.plot(data,data2)
        else:
            g.plot(data)
        #
        # Make a giffile
        #
        if not residue:
            residue='Pict1'
        psfile=residue+'.ps'
        if not os.path.isdir(dir):
            os.mkdir(dir)
 
        giffile=os.path.join(dir,residue+'.jpg')
        g.hardcopy(psfile, enhanced=1, color=1,mode='portrait')
        time.sleep(0.05)
        while 1:
            while not os.path.isfile(psfile):
                pass
            status=os.system('/usr/bin/convert "'+psfile+'" "'+giffile+'"')
            if os.path.isfile(giffile):
                break
            print 'Retry',
            #if status!=0:
            #    raise 'Non zero exit code from convert: ',status
        os.unlink(psfile)
        return giffile


    def makewwwpage(self,filename,giffiles,data,dpkas,title='Titration curves',table=1):
        #
        # Makes a www-page with all titration curves and a table of all dpka values
        #
        import os, string
        fd=open(filename,'w')
        fd.write('<HTML><HEAD><TITLE>'+title+'</TITLE></HEAD>\n')
        fd.write('<BODY BGCOLOR=EEEEEE text=black link=blue vlink="00007f">\n')
        fd.write('<h1>'+title+'</h1><BR><BR><HR>\n')
        if table==1:
            fd.write(string.join(self.insert_table(data,dpkas,giffiles)))
            fd.write('<BR><HR>')
        fd.write(string.join(self.insert_titration_curves(giffiles)))
        fd.write('</BODY></HTML>')
        fd.close()
        return

    def insert_table(self,data,dpka,giffiles):
        #
        # Makes a table of dpKa values and of the fit to the HH shape
        # data and dpka are dictionaries
        #
        lines=[]
        lines.append('Fit to Henderson-Hasselbalch equation. Slope should be close to 1\n')
        lines.append('<table bgcolor=lightblue cellspacing=2 border=5><tr><th colspan=2>Legend<tr><td><font color=pink>pink</font><td>Residue not present in the structure<br>\n')
        lines.append('<tr><Td><font color=yellow>yellow</font><td>Residue with undetermined pKa value\n')
        lines.append('<tr><td><font color=red>red</font><td>Residue with increased pKa value as compared to solution pKa value.\n')
        lines.append('<tr><td><font color=ligthgreen>green</font><td>Residue with decreased pKa value as compared to solution pKa value.</table>')
        lines.append('<TABLE BORDER=2 WIDTH=50% CELLSPACING=0 BGCOLOR=lightblue>\n')
        lines.append('<TH><center>Residue</center> <th><center>pKa</center><th><center> dpKa</center> <TH><center>Slope</center><TH><center>d</center> \n')
        residues=data.keys()
        residues.sort()
        for residue in residues:
            import string
            picfile=None
            for file in giffiles:
                if string.find(file,residue)!=-1:
                    picfile=file
            #
            # In case the data for this resiude does not exist
            #
            if not dpka[residue]['dpKa']:
                if dpka[residue]['pKa']:
                    dpka[residue]['dpKa']=0.0
                    
                else:
                    lines.append('<TR BGCOLOR="%s" align=center>' %'pink')
                    lines.append('<TD> %s <td> - <TD> - <TD> - <TD> -' %residue)
                    continue
            #
            # Do the colouring
            #
            nopka=None
            if abs(dpka[residue]['dpKa'])>50.0:
                #
                # This is if we were unable to determine the pKa value
                #
                colour='yellow'
                nopka=1
            elif abs(dpka[residue]['dpKa'])>1.0:
                colour='red'
            else:
                colour='lightgreen'
            lines.append('<TR BGCOLOR="%s" align=center>' %colour) 
            if picfile:
                link='<A HREF="%s" TARGET=BODY>%s</a>' %(picfile,residue)
            else:
                link='%s' %(residue)
            #
            # Make adjustments
            #
            if nopka:
                #
                # ok, so we couldn't determine the pKa value. Now lets find out if the residue is charged or uncharged
                #
                phs=self.tdata[residue].keys()
                neg=0
                pos=0
                for ph in phs:
                    if self.tdata[residue][ph]<-0.5:
                        neg=neg+1
                    elif self.tdata[residue][ph]>0.5:
                        pos=pos+1
                charge=0.0
                if neg!=0 and pos!=0:
                    raise "group is both negative and positive ",residue
                elif neg==0 and pos==0:
                    charge=0.0
                else:
                    if neg==0:
                        charge=1.0
                    else:
                        charge=-1.0
                #
                # Find out if the pKa is <0.0 or >20.0 ...
                #
                restype=residue[-3:]
                pka=None
                if restype=='ASP' or restype=='GLU' or restype=='TYR' or restype=='SER' or restype=='THR' or restype=='CYS':
                    if charge==0.0:
                        pka='> 20.0'
                    elif charge==-1.0:
                        pka='<  0.0'
                    else:
                        raise "Positive acid"
                else:
                    if charge==0.0:
                        pka='<  0.0'
                    elif charge==1.0:
                        pka='> 20.0'
                    else:
                        print restype,residue
                        raise "Negative base"
                lines.append('<TD> %s <td> %s <TD> %5.3f <TD> %5.3f <TD> %5.3f \n' %(link,pka,0.0,data[residue][0],data[residue][2]))
            else:
                lines.append('<TD> %s <td> %5.3f <TD> %5.3f <TD> %5.3f <TD> %5.3f\n' %(link,dpka[residue]['pKa'],dpka[residue]['dpKa'],data[residue][0],data[residue][2]))
        lines.append('</TABLE>\n\n')
        return lines
                     

    def insert_titration_curves(self,giffiles):
        #
        # insert all the titration curves
        #
        import os
        lines=[]
        count=0
        link=''
        picsperline=4
        for file in giffiles:
            if file[0]=="@":
                 continue
            linktext=os.path.split(file)[1][:-4]
            count=count+1
            if count<picsperline:
                lines.append('<A HREF="'+file+'"><IMG SRC="'+file+'" ALT='+linktext+' ALIGN=LEFT WIDTH=250 ></a>\n\n')
            else:
                lines.append('<A HREF="'+file+'"><IMG SRC="'+file+'" ALT='+linktext+' WIDTH=250 ></a><BR>\n\n')
            if count>1:
                link=link+'<SPACER SIZE=210><A HREF="'+file+'">'+linktext+'</a>\n'
            if count==1:
                link=link+'<spacer size=80><A HREF="'+file+'">'+linktext+'</a>\n'
            if count==picsperline:
                lines.append('<BR>'+link+'<BR><BR><HR NOSHADE><BR>\n\n')
                link=''
                count=0
        return lines
    

    def plot_stability(self,stabdata):
        #
        # This subroutine plots the contribution of each residue to
        # the stability of the protein. The routine uses gnuplot and is
        # only useful for very small proteins. For larger proteins the
        # plot is simply too complicated to analyse. In those cases I
        # suggest that you use the stability.dat and destability.dat
        # files created by analyse_stability in the pKanalyse class
        # wiht Excel or another plotting program.
        #
        import sys, Gnuplot, os, string, time
        stab={}
        destab={}
        prevres=''
        g=Gnuplot.Gnuplot()
        #
        # Add up all values so that we get a better plot
        #
        residues=stabdata.keys()
        residues.sort()
        for residue in residues:
            if not stab.has_key(residue):
                stab[residue]={}
            if not destab.has_key(residue):
                destab[residue]={}
            sys.stdout.flush()
            curve=stabdata[residue]
            for ph in curve.keys():
                if prevres!='':
                    if stabdata[residue][ph]<=0.0:
                        stab[residue][ph]=stab[prevres][ph]+stabdata[residue][ph]
                        destab[residue][ph]=destab[prevres][ph]
                    else:
                        stab[residue][ph]=stab[prevres][ph]
                        destab[residue][ph]=destab[prevres][ph]+stabdata[residue][ph]
                else:
                    if stabdata[residue][ph]<=0.0:
                        stab[residue][ph]=stabdata[residue][ph]
                        destab[residue][ph]=0.0
                    else:
                        stab[residue][ph]=0.0
                        destab[residue][ph]=stabdata[residue][ph]
            prevres=residue
        #
        # Now plot the data.
        #
        for stabdata in [stab,destab]:
            bigdata=[]
            max=-99.9
            min=99.9
            for residue in residues:
                plotdata=[]
                curve=stabdata[residue]
                phvals=curve.keys()
                phvals.sort()
                for ph in phvals:
                    plotdata.append(ph,curve[ph])
                    if curve[ph]>max:
                        max=curve[ph]
                    if curve[ph]<min:
                        min=curve[ph]        
                #bigdata.append(Gnuplot.Data(plotdata,title=residue,_with='lines'))
            #============================================================
            # Set min and max 
            #
            rmax=round(max)
            rmin=round(min)
            if rmax==rmin:
                if max>rmax:
                    rmax=rmax+1
                if min<rmin:
                    rmin=rmin-1
            if rmax==rmin:
                rmax=rmax+1
            rmax=max+2.0
            rmin=rmin-2.0
            #
            # Make the plot
            #
            g.title('Stability')
            g.xlabel('pH')
            g.ylabel('Free Energy')
            g('set yrange ['+str(rmin)+':'+str(rmax)+']')
            g.plotcmd = 'plot'
            g._clear_queue()
            for data in bigdata:
                g.itemlist.append(data)
            g.refresh()
            a=raw_input('Wait')
        return

            
        
#============================================================= 
            

class pKatools:

    def __init__(self):
        return

    def make_titration_curve(self,pKa,acid=1,residue='Hypothetical',pHstart=0.0,pHstop=20.0,pHstep=0.1):
        # Constructs a titration curve with the give pKa
        # Default is an acid constructed from pH 0.0 to 20.0 in 0.1 steps.
        dict={}
        for pH in range(100*pHstart,100*pHstop,100*pHstep):
            realpH=float(pH)/100.0
            charge=(1.0/(10.0**(realpH-pKa)+1.0))
            if acid==1:
                 charge=charge-1.0
            dict[realpH]=int(1000.0*charge)/1000.0
        dict['pKa']=pKa
        return {residue:dict}

    def make_titration_curve_2(self,pKa,acid,residue='hyp',pHvalues=[]):
        # Constructs a titration curve with the give pKa
        # Default is an acid constructed from pH 0.0 to 20.0 in 0.1 steps.
        dict={}
        import math
        for pH in pHvalues:
            charge=(1.0/(math.pow(10.0,(pH-pKa))+1.0))
            if acid==-1:
                 charge=charge-1.0
            else:
                charge=charge
            dict[pH]=charge
        dict['pKa']=pKa
        return {residue:dict}

    #
    # Get titratable residues in a PDB file
    #

    def get_titratable_residues(self,pdbfile):
        #
        # Call WHAT IF to get the titratable residues
        #
        import WI_tools, string, os
        pdbfile=os.path.join(os.getcwd(),pdbfile)
        titres=WI_tools.get_titratable_residues(pdbfile)
        group=[]
        for res in titres:
            res2=string.replace(res,'(','')
            res2=string.replace(res,')','')
            split=string.split(res2)
            #
            # Terminal residue?
            #
            term=None
            if split[-1][-4:]=='TERM':
                term=split[-1]
                split=split[:-1]
            if term:
                #
                # If a terminal we get rid of the residue name
                #
                resid=getWI_resid(res)
                resname=resid.split(':')[-1]
                resname=resid.replace(resname,'')
                group.append([resname+term,int(float(split[-1]))])
            else:
                group.append([getWI_resid(res),int(float(split[-1]))])
        return group


