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

import sys, os
import pKarun.pKa_utility_functions as pKaUT

class pKaIO:

    def __init__(self,rootfilename=None):
        self.newname=None
        self.backgr={}
        self.desolv={}
        self.matrix={}
        self.pka={}
        self.backgr_file=None
        self.desolv_file=None
        self.matrix_file=None
        self.pkafile=None
        self.titcurvfile=None
        #
        # Do we know any of the names?
        #
        if rootfilename:
            if rootfilename[-4:]!='.pka':
                self.backgr_file=rootfilename+'.BACKGR.DAT'
                self.desolv_file=rootfilename+'.DESOLV.DAT'
                self.matrix_file=rootfilename+'.MATRIX.DAT'
                self.pkafile=rootfilename+'.PKA.DAT'
                self.titcurvfile=rootfilename+'.TITCURV.DAT'
            else:
                self.backgr_file=rootfilename
                self.desolv_file=rootfilename
                self.matrix_file=rootfilename
                self.pkafile=rootfilename
                self.titcurvfile=rootfilename
        self.allfiles=[self.backgr_file,self.desolv_file,self.matrix_file,self.pkafile,self.titcurvfile]
        #
        # Do we have the files
        #
        self.assess_status()
        #
        # Define dictionaries that we will need all the time
        #
        import pKadata 
        self.acidbase=pKadata.acidbase
        self.modelpKas=pKadata.modelpKas
        #
        # We always use the new names
        #
        self.newname=1
        #
        # done
        #
        return
        

    #
    # --------------------
    #

    def assess_status(self):
        #
        # Do we have a completed pKa calculation for the PDB file?
        #
        import os
        self.file_status={}
        self.calculation_completed=1
        for file in self.allfiles:
            if not file:
                self.calculation_completed=None
                self.file_status[file]=None
                continue
            if os.path.isfile(file):
                self.file_status[file]=1
            else:
                self.calculation_completed=None
                self.file_status[file]=None
        return self.calculation_completed

    #
    # -------------------------------
    #
    
    def read_MEADfile(self,lines,mode='pKa'):
        """Parse a MEAD pKa calculation file"""
        # Some residues have different names
        MEADtranslate={'ASH1':'ASP','GLH1':'GLU'}
        #
        count=0
        while count<len(lines):
            line=lines[count]
            if line.strip()=='<****CALCULATION RESULTS>' and mode=='pKa':
                pKas={}
                while line[:10]!='--------------'[:10] or count==len(lines):
                    count=count+1
                    line=lines[count]
                if line[:10]!='---------------'[:10]:
                    raise Exception('Could not parse MEAD file')
                #
                count=count+1
                line=lines[count]
                while line.strip()!='':
                    sp=line.split()
                    residue=sp[0].split('-')
                    import string
                    residue='%s:%s:%s' %('',string.zfill(int(residue[1]),4),residue[0])
                    modelpk=float(sp[1])
                    desolv=float(sp[3])
                    backgr=float(sp[2])
                    intpka=float(sp[4])
                    pka=float(sp[5])
                    site=pka-intpka
                    titgroup=residue.split(':')[2]
                    if MEADtranslate.has_key(titgroup):
                        titgroup=MEADtranslate[titgroup]
                    pKas[residue]={'pKa':pka,'modelpK':modelpk,'desolv':desolv,'backgr':backgr,'delec':site,
                          'acidbase':self.acidbase[titgroup],'intpka':intpka}
                    count=count+1
                    line=lines[count]
                return pKas
            #
            # Code for reading the titration curves
            #
            elif line.strip()=='<****SITE_INFOS>' and mode=='titcurv':
                self.titdata={}
                count=count+1
                titcurvs=eval(lines[count])
                for item in titcurvs:
                    residue=item['site'][0].split('-')
                    import string
                    residue='%s:%s:%s' %('',string.zfill(int(residue[1]),4),residue[0])
                    tc=item['tit_curve']
                    this_tc={}
                    #
                    import pKadata
                    acidbase=pKadata.acidbase
                    for pH,charge in tc:
                        titgroup=residue.split(':')[2]
                        if MEADtranslate.has_key(titgroup):
                            titgroup=MEADtranslate[titgroup]
                        if not acidbase.has_key(titgroup):
                            raise Exception('Unknown titratable group: %s\nPlease add this titratable group to pKaTool/pKaData.py' %titgroup)
                        this_tc[float(pH)]=float(acidbase[titgroup])*float(charge)
                    self.titdata[residue]=this_tc.copy()
                return self.titdata
            #
            # Code for reading the interaction matrix
            #
            elif line.strip()=='<****SITE_INFOS>' and mode=='matrix':
                self.matrix={}
                count=count+1
                infodict=eval(lines[count])
                #
                # Get all the titgroups first
                #
                titgroups=[]
                for item in infodict:
                    residue=item['site'][0].split('-')
                    import string
                    residue='%s:%s:%s' %('',string.zfill(int(residue[1]),4),residue[0])
                    titgroups.append(residue)
                #
                # Now we can read teh matrix
                #
                for item in infodict:
                    residue=item['site'][0].split('-')
                    import string
                    residue='%s:%s:%s' %('',string.zfill(int(residue[1]),4),residue[0])
                    self.matrix[residue]={}
                    ints=item['interactions']
                    count=0
                    for interaction in ints:
                        import pKadata
                        acidbase=pKadata.acidbase
                        titgroup1=residue.split(':')[2]
                        if MEADtranslate.has_key(titgroup1):
                            titgroup1=MEADtranslate[titgroup1]
                        group1=acidbase[titgroup1]
                        #
                        titgroup2=titgroups[count].split(':')[2]
                        if MEADtranslate.has_key(titgroup2):
                            titgroup2=MEADtranslate[titgroup2]
                        group2=acidbase[titgroup2]
                        #
                        self.matrix[residue][titgroups[count]]=[331.842/0.592*interaction*float(group1)*float(group2),0.0,0.0,0.0] # Multiply to get the sign right
                        # 331.842 is a conversion factor to go from internal MEAD energies to kcal/mol
                        count=count+1
                return self.matrix
            #
            # Read the coordinates from the file
            #
            elif line.strip()=='<****ATOMS>' and mode=='pdb':
                count=count+2
                line=lines[count]
                pdblines=[]
                while line.strip()!='' and line.strip()!='<****CONNECT>':
                    sp=line.split()
                    number=int(sp[0])
                    resnumber=sp[1]
                    resname=sp[2]
                    #chainid=sp[3]
                    chainid=' '
                    atomname=sp[4]
                    X=float(sp[6])
                    Y=float(sp[7])
                    Z=float(sp[8])
                    pdbline='ATOM  %5i %4s %-4s%1s%4s    %8.3f%8.3f%8.3f %5s %5s' %(number,atomname,resname,chainid,resnumber,X,Y,Z,0,0)
                    pdblines.append(pdbline)
                    count=count+1
                    line=lines[count]
                return pdblines
            count=count+1
        raise Exception('Could not parse MEAD file')
        
    #
    # --------------
    #
    
    def readpka(self,filename=None):
        """
        # This function reads a WHAT IF pKa file and creates a dictionary with the
        # following format: self.pka={<residue1>:{'pKa':<pKa>,'modelpK':<model pKa value>,
        # 'desolv':<dpK due to desolvation>,'backgr':<dpK due to background int.>,
        # 'delec':<dpK due to site-site interactions>},<residue2>:.....}
        """
        if not filename:
            filename=self.pkafile
        import os, string
        if not filename:
            print filename,'is not a filename'
            os._exit(0)
        if not os.path.isfile(filename):
            raise Exception('File does not exist: %s' %filename)
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        #
        # Parse
        #
        if string.strip(lines[0]).lower()=='WHAT IF pKa File'.lower():
            format='WHAT IF'
        elif string.lower(string.strip(lines[0]))==string.lower('pdb2pka pKa File'):
            format='WHAT IF'
        elif lines[2].strip()=='#pKa initialization':
            return self.read_MEADfile(lines,mode='pKa')
        else:
            raise Exception('Unknown format')
        if string.lower(string.strip(lines[1]))!=string.lower('Format 1.0'):
            raise Exception('unknown format: %s' %str(lines[1]))
        # Next line is text
        linenumber=3
        pKa={}
        done=None
        #
        # Define counter for keeping track of terminals
        #
        Nterm=0
        while not done:
            residue=string.split(lines[linenumber])
            if format=='WHAT IF':
                nline=''
                if string.lower(residue[0])==string.lower('TERMINAL'):
                    Nterm=abs(Nterm-1)
                    nline=string.replace(lines[linenumber+1],'(','')
                else:
                    nline=string.replace(lines[linenumber],'(','')
                nline=string.replace(nline,')','')
                list=string.split(nline)
                if len(list[3])==1:
                    resid=string.zfill(list[2],4)+string.upper(list[1])+list[3]
                else:
                    resid=string.zfill(list[2],4)+string.upper(list[1])
                    chainid=None
                if string.lower(residue[0])==string.lower('TERMINAL'):
                    linenumber=linenumber+1
                    residue=string.split(lines[linenumber])
                    if len(residue[5])==1 and residue[5] in string.letters:
                        pka=string.atof(residue[6])
                    else:
                        pka=string.atof(residue[5])
                    if len(residue)>6:
                        if len(residue[5])==1 and residue[5] in string.letters:
                            modelpk=string.atof(residue[7])
                            dpk=string.atof(residue[8])
                            desolv=string.atof(residue[0])
                            backgr=string.atof(residue[10])
                            site=string.atof(residue[11])
                        else:
                            index=6
                            if len(residue[2])>1:
                                index=index-1
                            modelpk=string.atof(residue[index])
                            dpk=string.atof(residue[index+1])
                            desolv=string.atof(residue[index+2])
                            backgr=string.atof(residue[index+3])
                            site=string.atof(residue[index+4])
                    else:
                        modelpk='Not in file'
                        dpk='Not in file'
                        desolv='Not in file'
                        backgr='Not in file'
                        site='Not in file'

                    #residue='T'+string.zfill(residue[0],4)+residue[1]
                    residue='T'+resid
                else:
                    if len(residue[5])==1 and residue[5] in string.letters:
                        pka=string.atof(residue[6])
                    else:
                        pka=string.atof(residue[5])
                    #pka=string.atof(residue[5])
                    if len(residue)>6:
                        if len(residue[5])==1 and residue[5] in string.letters:
                            modelpk=string.atof(residue[7])
                            dpk=string.atof(residue[8])
                            desolv=string.atof(residue[9])
                            backgr=string.atof(residue[10])
                            site=string.atof(residue[11])
                        else:
                            index=6
                            if len(residue[2])>1:
                                index=index-1
                            modelpk=string.atof(residue[index])
                            dpk=string.atof(residue[index+1])
                            desolv=string.atof(residue[index+2])
                            backgr=string.atof(residue[index+3])
                            site=string.atof(residue[index+4])
                    else:
                        modelpk='Not in file'
                        dpk='Not in file'
                        desolv='Not in file'
                        backgr='Not in file'
                        site='Not in file'
                    #residue=string.zfill(residue[0],4)+residue[1]
                    residue=resid
            elif format=='pdb2pka':
                pka=string.atof(residue[1])
                modelpk=string.atof(residue[2])
                dpk=string.atof(residue[3])
                desolv=string.atof(residue[4])
                backgr=string.atof(residue[5])
                site=string.atof(residue[6])
                residue=residue[0]
            #
            # Reformat the residue names if we are asked for the new name
            # structure
            #
            if self.newname:
                residue=pKaUT.reformat_name(residue,Nterm,format)
            #
            # Make sure that non-determined pKa values are marked appropriately
            #
            if pka==0.0:
                pka=None
            #
            # Construct dictionary
            #
            pKa[residue]={'pKa':pka,'modelpK':modelpk,'desolv':desolv,'backgr':backgr,'delec':site,
                          'acidbase':self.acidbase[residue.split(':')[2]],'intpka':modelpk+desolv+backgr}
            linenumber=linenumber+1
            if string.strip(string.lower(lines[linenumber]))==string.lower('End of file'):
                done=1
        #
        # Change old style terminal names 'A:0010:ASP:CTERM' to new style 'A:0010:CTERM'
        #
        import string
        for group in pKa.keys():
            sp=group.split(':')
            if sp[-1]=='CTERM' or sp[-1]=='NTERM':
                if sp[-2][0] in string.letters:
                    pKa['%s:%s:%s' %(sp[0],sp[1],sp[-1])]=pKa[group]
                    del pKa[group]
        #
        # Return data
        #
        self.pka=pKa
        return self.pka

    #
    # ----
    #




    #
    # -----
    #

    def write_pka(self,filename,data=None,format='WHAT IF'):
        """Write a PKA.DAT file containing all info on the calculations"""
        #
        # Get the data
        #
        if not data:
            data=self.pka
        #
        # Write the file
        #
        fd=open(filename,'w')
        fd.write('%s pKa File\n' %format)
        fd.write('Format 1.0\n')
        fd.write('      Group            pKa value  Model pK    dpK     dDesolv    dBack   dElec\n')
        groups=data.keys()
        groups.sort()
        written={}
        #
        # ---------
        #
        for group in groups:
            if data.has_key(group):
                this_data=data[group]
                fd.write('%15s      %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  \n' %(self.WI_res_text(group,format),
                                                                                    this_data['pKa'],
                                                                                    this_data['modelpK'],
                                                                                    this_data['pKa']-this_data['modelpK'],
                                                                                    this_data['desolv'],
                                                                                    this_data['backgr'],
                                                                                    this_data['delec']))
                written[group]=1
        fd.write('End of file\n')
        fd.close()
        return
                         
    #
    # -------------------------------
    #

    def read_titration_curve(self,filename=None):
        return self.readtitcurv(filename)

    def readtitcurv(self,filename=None):
        #
        # Syntax: readtitcurv(self,<titration curve filename>)
        # This function reads a WHAT IF titration curve file and
        # creates self.titdata, which is a dictionary:
        # self.titdata={<residue>:{'pKa':<pka>,<ph1>:<charge1>,<ph2>:<charge2>.....}}
        #
        if not filename:
            filename=self.titcurvfile
        import os, string
        if not os.path.isfile(filename):
            raise Exception('File does not exist: %s' %filename)
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        #
        # Parse
        #
        if string.lower(string.strip(lines[0]))==string.lower('WHAT IF Titration Curve File'):
            format='WHAT IF'
        elif string.lower(string.strip(lines[0]))==string.lower('pdb2pka Titration Curve File'):
            format='WHAT IF'
        elif lines[2].strip()=='#pKa initialization':
            return self.read_MEADfile(lines,mode='titcurv')
        else:
            raise Exception('Unknown format')
        if string.lower(string.strip(lines[1]))!=string.lower('Format 1.0'):
            raise Exception('unknown format: '+str(lines[1]))
        #
        # Get the pH start, stop and step
        #
        phvals=string.split(lines[2])
        phstart=string.atof(phvals[0])
        phend=string.atof(phvals[1])
        phstep=string.atof(phvals[2])
        titdata={}
        linenumber=3
        done=0
        terms=['NTERM','CTERM']
        term_count=-1
        while not done:
            Term=False
            if string.strip(lines[linenumber])=='TERMINAL GROUP:':
                linenumber=linenumber+1
                Term=1
            residue=pKaUT.getWI_resid2(lines[linenumber],format)
            if Term:
                term_count=term_count+1
                #
                # Get rid of the residue name and just add the term
                #
                aaname=residue.split(':')[-1]
                residue=residue.replace(aaname,terms[term_count])
                if term_count==1:
                    term_count=-1
            pKa=float(string.split(lines[linenumber])[-1])
            linenumber=linenumber+1
            #
            # --------------
            #
            charge={'pKa':pKa}
            for pH in range(int(100*phstart),int(100*phend+100*phstep),int(100*phstep)):
                rpH=float(pH)/100.0
                line=string.split(lines[linenumber])
                if string.atof(line[0])==rpH:
                    charge[rpH]=string.atof(line[1])
                    linenumber=linenumber+1
            titdata[residue]=charge
            linenumber=linenumber+1
            if string.strip(string.lower(lines[linenumber]))==string.lower('End of file'):
                done=1
        self.titdata=titdata
        return self.titdata

    #
    # ----------------------
    #

    def write_titration_curve(self,filename,data,format='WHAT IF'):
        #
        # This function creates a WHAT IF titration curve file.
        # data is a dictionary:
        # data={<residue>:{'pKa':<pka>,<ph1>:<charge1>,<ph2>:<charge2>.....}}
        #
        # Extract some data from the dictionary
        #
        residues=data.keys()
        phvals=data[residues[0]].keys()
        phvals.sort()
        for residue in residues:
            newpHvals=data[residue].keys()
            newpHvals.sort()
            if newpHvals!=phvals:
                print phvals
                print newpHvals
                raise Exception('Dictionary does not contain identical pH values')
        print 'Done'
        #
        # Check that a pKa value is in the pH values
        #
        for residue in residues:
            if not data[residue].has_key('pKa'):
                data[residue]['pKa']=0.0
        #
        # Find the pH-start, stop and step
        #
        phvals=data[residues[0]].keys()
        phvals.sort()
        phstart=phvals[0]
        phstop=phvals[-2]
        phstep=phvals[1]-phstart
        lastval=phstart
        #
        # Write the file
        #
        import os, string
        fd=open(filename,'w')
        #
        # Write header
        #
        fd.write('%s Titration Curve File\n' %format)
        fd.write('Format 1.0\n')
        #
        # Start pH, end pH, pH step
        #
        fd.write('%6.3f %7.3f %6.3f\n' %(phstart,phstop,phstep))
        residues=data.keys()
        residues.sort()
        for residue in residues:
            fd.write('%s      %7.4f\n' %(self.WI_res_text(residue,format),float(data[residue]['pKa'])))
            for ph in phvals:
                if ph=='pKa':
                    continue
                fd.write('%.2f  %.3f\n' %(float(ph),float(data[residue][ph])))
            fd.write('------------------------------------------\n')
        fd.write('End of file\n')
        #
        # Close file
        #
        fd.close()
        #
        # This will never work without a template file.
        # It is not worth the trouble to reconstruct WHAT IFs residue identifier line
        #
        return

    #
    # ----------------------------------
    #

    def read_matrix(self,filename=None):
        """
        # This subroutine read a MATRIX file
        """
        if not filename:
            if self.matrix_file:
                filename=self.matrix_file
            else:
                raise Exception('No matrix filename given')
        #
        import os, string
        if not os.path.isfile(filename):
            raise Exception("File not found: "+filename)
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        #
        # Initialise dictionary
        #
        self.matrix={}
        #
        # Read title lines
        #
        if string.lower(string.strip(lines[0]))==string.lower('WHAT IF Interaction Matrix File'):
            format='WHAT IF'
        elif string.lower(string.strip(lines[0]))==string.lower('pdb2pka Interaction Matrix File'):
            format='WHAT IF'
        elif lines[2].strip()=='#pKa initialization':
            return self.read_MEADfile(lines,mode='matrix')
        else:
            raise Exception('Unknown format')
        if not string.strip(lines[1])=='Format 1.0':
            raise Exception('Wrong format'+str(lines[1]))
        x=1
        done=None
        partners=None
        Nterm=0
        while not done:
            x=x+1
            #
            # Read first line for this partner residue
            #
            if format=='WHAT IF':
                term=None
                if string.strip(lines[x])=='TERMINAL GROUP:':
                    term=1
                    Nterm=abs(Nterm-1)
                    x=x+1
                nline=string.replace(lines[x],'(','')
                nline=string.replace(nline,')','')
                resid=pKaUT.getWI_resid2(lines[x])
                if term:
                    #
                    # Replace the residue name by the terminal name
                    #
                    resname=resid.split(':')[-1]
                    d={1:'NTERM',0:'CTERM'}
                    resid=resid.replace(resname,d[Nterm])
                #
                # Get the number of interaction partners
                #
                np=int(float(nline.split()[-1]))
                if not partners:
                    partners=np
                else:
                    if partners!=np:
                        raise Exception('Number of partners changes: %d' %np)
                self.matrix[resid]={}
                #
                # Now read all the interactions with the partners
                #
                Nterm_partner=0
                for count in range(partners):
                    x=x+1
                    term2=None
                    if string.strip(lines[x])=='TERMINAL GROUP:':
                        Nterm_partner=abs(Nterm_partner-1)
                        term2=1
                        x=x+1
                    nline=string.replace(lines[x],'(','')
                    nline=string.replace(nline,')','')
                    list=string.split(nline)
                    partid=pKaUT.getWI_resid2(lines[x])
                    if term2:
                        #
                        # Replace the residue name by the terminal name
                        #
                        resname=partid.split(':')[-1]
                        d={1:'NTERM',0:'CTERM'}
                        partid=partid.replace(resname,d[Nterm_partner])
                    #
                    # New name?
                    #
                    chacha=string.atof(string.strip(list[-1]))
                    x=x+1
                    i2=float(lines[x])
                    x=x+1
                    i3=float(lines[x])
                    x=x+1
                    i4=float(lines[x])
                    energies=[chacha,i2,i3,i4]
                    self.matrix[resid][partid]=energies
                    term2=None
            elif format=='pdb2pka':
                #
                # pseudo-pdb2pka format
                #
                list=lines[x].split()
                resid=string.strip(list[0])
                self.matrix[resid]={}
                partners=int(float(list[-1]))
                #
                # Now read all the interactions with the partners
                #
                for count in range(partners):
                    x=x+1
                    list=lines[x].split()
                    partid=list[0]
                    chacha=string.atof(string.strip(list[-1]))
                    x=x+1
                    i2=string.atof(lines[x])
                    x=x+1
                    i3=string.atof(lines[x])
                    x=x+1
                    i4=string.atof(lines[x])
                    energies=[chacha,i2,i3,i4]
                    self.matrix[resid][partid]=energies
            x=x+1
            if string.strip(lines[x+1])=='End of file':
                done=1
        return self.matrix

    #
    # -------------------------
    #

    def write_matrix(self,filename,format='WHAT IF'):
        #
        # Writes an interaction energy matrix
        #
        fd=open(filename,'w')
        fd.write('%s Interaction Matrix File\n' %format)
        fd.write('Format 1.0\n')
        groups=self.matrix.keys()
        groups.sort()
        num_groups=len(groups)
        count=0
        written={}
        newgroups=[]
        for group in groups:
            if group[0]=='T':
                newgroups.append(group[1:])
            else:
                newgroups.append(group)
        newgroups.sort()
        self.newgroups=newgroups[:]
        #
        # ---------
        #
        for group in newgroups:
            if self.matrix.has_key(group):
                fd.write('%s      %7.4f\n' %(self.WI_res_text(group,format),float(num_groups)))
                self.write_section(group,fd,format)
                written[group]=1
                #
                # Is there a terminal group associated with this residue?
                #
                if self.matrix.has_key('T'+group):
                    fd.write('%s      %7.4f\n' %(self.WI_res_text('T'+group),float(num_groups)))
                    self.write_section('T'+group,fd,format)
                    written['T'+group]=1
            else:
                if self.matrix.has_key('T'+group) and not written.has_key('T'+group):
                    fd.write('%s      %7.4f\n' %(self.WI_res_text('T'+group,format),float(num_groups)))
                    self.write_section('T'+group,fd,format)
                    written['T'+group]=1
        fd.write('End of file\n')
        return

    #
    # -----
    #

    def write_pdb2pka_matrix(self,filename,matrix):
        """Write an interaction energy matrix in pdb2pka format"""
        """At the moment, we just reformat and write a WHAT IF file"""
        self.matrix={}
        for group1 in matrix.keys():
            self.matrix[group1.uniqueid]={}
            for tit1 in matrix[group1].keys():
                for state1 in matrix[group1][tit1].keys():
                    sub_m=matrix[group1][tit1][state1]
                    for group2 in sub_m.keys():
                        if not self.matrix[group1.uniqueid].has_key(group2.uniqueid):
                            self.matrix[group1.uniqueid][group2.uniqueid]=[]
                        for tit2 in sub_m[group2].keys():
                            for state2 in sub_m[group2][tit2].keys():
                                self.matrix[group1.uniqueid][group2.uniqueid].append(sub_m[group2][tit2][state2])
        for group1 in self.matrix.keys():
            for group2 in self.matrix[group1].keys():
                sum=0.0
                for val in self.matrix[group1][group2]:
                    sum=sum+val
                self.matrix[group1][group2]=[sum,0.0,0.0,0.0]
        self.write_matrix(filename,format='pdb2pka')
        return
   
    #
    # ------------------------
    #

    def write_section(self,group,fd,format):
        groups_tmp=self.matrix[group].keys()
        groups2=[]
        for group_x in groups_tmp:
            if group_x[0]=='T':
                groups2.append(group_x[1:])
            else:
                groups2.append(group_x)
        groups2.sort()
        written={}
        for group2 in groups2:
            if self.matrix[group].has_key(group2):
                fd.write('%s      %7.4f\n' %(self.WI_res_text(group2,format),self.matrix[group][group2][0]))
                fd.write('%7.4f\n%7.4f\n%7.4f\n'%(self.matrix[group][group2][1],self.matrix[group][group2][2],self.matrix[group][group2][3]))
                written[group2]=1
                #
                # Is there a terminal group associated with this residue?
                #
                if self.matrix[group].has_key('T'+group2):
                    fd.write('%s      %7.4f\n' %(self.WI_res_text('T'+group2,format),self.matrix[group]['T'+group2][0]))
                    fd.write('%7.4f\n%7.4f\n%7.4f\n'%(self.matrix[group]['T'+group2][1],self.matrix[group]['T'+group2][2],self.matrix[group]['T'+group2][3]))
                    written['T'+group2]=1
            else:
                if self.matrix[group].has_key('T'+group2) and not written.has_key('T'+group2):
                    fd.write('%s      %7.4f\n' %(self.WI_res_text('T'+group2,format),self.matrix[group]['T'+group2][0]))
                    fd.write('%7.4f\n%7.4f\n%7.4f\n'%(self.matrix[group]['T'+group2][1],self.matrix[group]['T'+group2][2],self.matrix[group]['T'+group2][3]))
                    written['T'+group2]=1 
        fd.write('--------------------------------------------\n')
        return

    #
    # ------------------------------
    #

    def read_desolv(self,filename=None):
        if not filename:
            if self.desolv_file:
                filename=self.desolv_file
            else:
                raise Exception('No desolv filename given')
        #
        #
        # This subroutine reads a DESOLV file
        #
        import os, string
        if not os.path.isfile(filename):
            raise Exception("File not found:%s" %filename)
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        #
        # Initialise dictionary
        #
        self.desolv={}
        #
        # Read title lines
        #
        if string.strip(lines[0])=='WHAT IF Desolvation Energy File' and string.strip(lines[1])=='Format 1.0':
            format='WHAT IF'
        elif string.strip(lines[0])=='pdb2pka Desolvation Energy File':
            format='pdb2pka'
        else:
            raise Exception('Unknown format:'+string.strip(lines[0]))
        #
        # Call the generic read routine
        #
        self.read_WIfile(lines,self.desolv,format)
        return self.desolv
    
    #
    # -----------------------------
    #

    def read_backgr(self,filename=None):
        #
        # This subroutine reads a BACKGR file
        #
        if not filename:
            if self.backgr_file:
                filename=self.backgr_file
            else:
                raise Exception('No matrix filename given')
        #
        import os, string
        if not os.path.isfile(filename):
            raise Exception("File not found: %s" %filename)
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        #
        # Initialise dictionary
        #
        self.backgr={}
        #
        # Read title lines
        #
        if string.strip(lines[0])=='WHAT IF Background Energy File':
            format='WHAT IF'
        elif string.strip(lines[0])=='pdb2pka Background Energy File':
            format='pdb2pka'
        else:
            raise Exception('Unknown format:'+string.strip(lines[0]))
        #
        # Call the generic read routine
        #
        self.read_WIfile(lines,self.backgr,format)
        return self.backgr

    #
    # ----------------------------------
    #


    def read_WIfile(self,lines,dict,format):
        """
        # Read a DESOLV file or a BACKGR file
        """
        import string
        x=1
        done=None
        Nterm=0
        while not done:
            x=x+1
            #
            # Read first line for this residue
            #
            term=None
            if lines[x].strip()=='TERMINAL GROUP:':
                Nterm=abs(Nterm-1)
                term=True
                x=x+1
            nline=string.replace(lines[x],'(','')
            nline=string.replace(nline,')','')
            list=string.split(nline)
            resid=pKaUT.getWI_resid2(lines[x])
            if term:
                #
                # Replace the residue name by the terminal name
                #
                resname=resid.split(':')[-1]
                d={1:'NTERM',0:'CTERM'}
                resid=resid.replace(resname,d[Nterm])
            #
            dict[resid]=string.atof(list[-1])
            if string.strip(lines[x+1])=='End of file':
                done=1
        return dict

    #
    # ----------------------
    #

    def write_desolv(self,filename,format='WHAT IF'):
        #
        # Writes the desolvation file
        #
        fd=open(filename,'w')
        fd.write('%s Desolvation Energy File\n' %format)
        fd.write('Format 1.0\n')
        groups=self.desolv.keys()
        groups.sort()
        written={}
        #
        # ---------
        #
        for group in groups:
            if self.desolv.has_key(group):
                fd.write('%s      %7.4f\n' %(self.WI_res_text(group,format),float(self.desolv[group])))
                written[group]=1
        fd.write('End of file\n')
        return

    #
    # ----------------------
    #

    def write_backgr(self,filename,format='WHAT IF'):
        #
        # Writes the background interaction energy file
        #
        fd=open(filename,'w')
        fd.write('%s Background Energy File\n' %format)
        fd.write('Format 1.0\n')
        groups=self.backgr.keys()
        groups.sort()
        written={}
        #
        # ---------
        #
        for group in groups:
            if self.backgr.has_key(group):
                fd.write('%s      %7.4f\n' %(self.WI_res_text(group,format),float(self.backgr[group])))
                written[group]=1
        fd.write('End of file\n')
        return

    #
    # ----------------------
    #
    
    def write_pKacoop_submission(self,filename='submit.pkac',prediction_name='TEST',groups=[]):
        """Write a submission file for the pKa cooperative"""
        if groups==[]:
            return
        groups.sort()
        tcs=self.readtitcurv()
        fd=open(filename,'w')
        for group in groups:
            fd.write('PREDICTION: %s\n' %prediction_name)
            number=group.split(':')
            name='%s %s %s' %(number[2],number[0],number[1])
            fd.write('TITGROUP: %s\n' %name)
            pHs=tcs[group].keys()
            pHs.sort()
            pHs.remove('pKa')
            for pH in pHs:
                fd.write('%5.2f, %5.2f\n' %(pH,tcs[group][pH]))
            fd.write('END\n')
        fd.close()
            
    #
    # ----------------------
    #

    def WI_res_text(self,residue,format):
        """Constructs the WHAT IF residue ID line
        f.ex. 1 LYS  (   1  ) from 0001LYS.
        Function works with new names."""
        if format=='WHAT IF':
            if not self.newname:
                # Old names
                import string
                terminal=None
                if residue[0]=="T":
                    terminal=1
                    residue=residue[1:]
                number=string.atoi(residue[:4])
                residue=residue[4:]
                if len(residue)>3:
                    chainid=residue[-1]
                    residue=residue[:-1]
                else:
                    chainid=' '
                line='%4d %3s  (%4d  ) %1s' %(number,residue,number,chainid)
                if terminal:
                    line='TERMINAL GROUP:\n'+line
                return line
            else:
                #
                # New names
                #
                terminal=None
                split=residue.split(':')
                if split[-1]=="CTERM" or split[-1]=='NTERM':
                    terminal=1
                try:
                    number=int(split[1])
                    residue=split[2]
                    chainid=split[0]
                    line='%4d %3s  (%4d  ) %1s' %(number,residue,number,chainid)
                    if terminal:
                        line='TERMINAL GROUP:\n'+line
                    return line
                except:
                    return residue
        elif format=='pdb2pka':
            #
            # PDB2PKA format
            #
            terminal=None
            split=residue.split(':')
            if split[1]=='NTR' or split[1]=='CTR':
                terminal=1
            #
            res=split[0]
            res=res.split('_')
            residue=res[0]
            chainid=res[1]
            #
            # For Design_pKa
            #
            chainid=''
            #
            number=int(res[2])
            line='%4d %3s  (%4d  ) %1s' %(number,residue,number,chainid)
            if terminal:
                line='TERMINAL GROUP:\n'+line
            return line
        else:
            raise Exception('Unknown format:'+format)


