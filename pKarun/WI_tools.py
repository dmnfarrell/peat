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

import WHAT_IF, string, os, pKa_ostools


def get_sequence(pdbfile):
    import os,string
    command='getmol '+os.path.split(pdbfile)[1]+' \n \n set \n %soupir \n tot \n seq.pir \n '
    logfile,files=RunWI(command,['seq.pir'],None,[pdbfile])
    try:
        return files['seq.pir']
    except:
        print 'Could not get sequence for ',pdbfile
        raise 'Error'

def superimpose(mol1,mol2):
    #
    # Superimpose mol2 on mol1
    # everything is superimposed
    #
    mol1s=os.path.split(mol1)[1]
    mol2s=os.path.split(mol2)[1]
    alignedfile=mol2s[:-3]+'aln.pdb'
    #Prepare commandline for WHATIF.
    com='getmol '+mol1s+'\n mol1 \n getmol '+mol2s+' mol2 \n %PASTML'
    com=com+' suppos \n params\n rmserr 2 maxerr 4 '
    com=com+' exit\n motif\n Y\n m1\n m2\n N\n apply mol2 0 \n'
    com=com+' soup makmol \n \n '+alignedfile+'\n '
    com=com+'No comments \n \n mol2\n 0\n'
    #here WHAT IF is called.
    output,files=RunWI(com,[alignedfile],None,[mol1,mol2])
    return output,files

def center_molecule(pdbfile):
    centerfile=pdbfile[:-4]+'cent.pdb'
    com='getmol '+pdbfile+' \n \n minmax tot 0 \n Y \n'
    com=com+' %makmol '+pdbfile+' \n '+centerfile+' \n \n tot 0 \n'
    output=RunWI(com)
    return centerfile

def calculate_phimap(pdbfile,phifile,parameters=None,deletewater=1,crgfile=None,radfile=None):
    home=os.getcwd()
    #
    # Make temporary directory
    #
    tempfile.tempdir=home
    dir=tempfile.mkdtemp()
    #
    # Make symlink to TOPOLOGY.H
    #
    os.link(pdbfile,os.path.join(dir,pdbfile))
    os.symlink(os.path.join(WHAT_IF.dbdata,'TOPOLOGY.H'),os.path.join(dir,'TOPOLOGY.H'))
    #
    # Make links to the charge file and the radius file if nonstandard radfiles are
    # requested.
    #
    if crgfile:
        if os.path.isfile(os.path.join(WHAT_IF.dbdata,crgfile)):
            os.symlink(os.path.join(WHAT_IF.dbdata,crgfile),os.path.join(dir,'DELCRG.DAT'))
        else:
            raise 'Charge file does not exist in '+WHAT_IF.dbdata
    if radfile:
        if os.path.isfile(os.path.join(WHAT_IF.dbdata,radfile)):
            os.symlink(os.path.join(WHAT_IF.dbdata,radfile),os.path.join(dir,'DELRAD.DAT'))
        else:
            raise 'Radius file does not exist in '+WHAT_IF.dbdata
    #
    # cd to dir and run WHAT IF
    #
    os.chdir(dir)
    command='getmol '+pdbfile+' \n \n '
    if deletewater==1:
        command=command+' %DELWAT \n '
    command=command+' %setwif 339 1 \n %addhyd tot 0 \n tot 0\n  '
    #
    # if any parameters were given the set them
    #
    if parameters:
        for key in parameters.keys():
            command=command+' setwif '+str(key)+' '+str(parameters[key])+' \n' 

    #
    # complete the command
    #
    command=command+' %setscg \n %rundel tot 0\n '+phifile+' \n '
    output=RunWI(command)
    os.chdir(home)
    os.link(os.path.join(dir,phifile),phifile)
    pKa_ostools.rmr(dir)
    return


#
# ----
#

def calculate_dipole(pdbfile):
    if os.path.isfile('TOPOLOGY.H'):
        os.unlink('TOPOLOGY.H')
    if os.path.isfile('DELRAD.DAT'):
        os.unlink('DELRAD.DAT')
    if os.path.isfile('DELCRG.DAT'):
        os.unlink('DELCRG.DAT')
    os.symlink(os.path.join(WHAT_IF.dbdata,'TOPOLOGY.H'),'TOPOLOGY.H')
    os.symlink(os.path.join(WHAT_IF.dbdata,'OPLS.CRG_noligand'),'DELCRG.DAT')
    os.symlink(os.path.join(WHAT_IF.dbdata,'OPLS.RAD'),'DELRAD.DAT')
    com=' getmol '+pdbfile+' \n Y \n \n %delwat \n %corall N \n \nsetwif 339 1 \n %addhyd tot 0 \n tot 0 \n'
    com=com+' setwif 1828 1 \n setwif 1829 0 \n setwif 1830 150 \n %setscg \n %eletot \n'
    output=RunWI(com)
    os.unlink('TOPOLOGY.H')
    os.unlink('DELRAD.DAT')
    os.unlink('DELCRG.DAT')
    return output

def do_fulchk(pdbfile):
    """Run a fulchk on a pdbfile"""
    import tempfile, os
    topdir=os.getcwd()
    tdir=tempfile.mkdtemp(dir=topdir)
    import shutil
    pdbname=os.path.join(tdir,os.path.split(pdbfile)[1])
    shutil.copyfile(pdbfile,pdbname)
    #
    os.chdir(tdir)
    newpdb=pdbname+'_delhyd'
    # First delete hydrogens
    com =' $pwd \n $ls \n getmol %s \n set1 \n 1 \n Y \n 1 \n %%dellig \n soup \n %%makmol \n %s \n %s \n \n tot 0 \n' %(pdbname,pdbname,newpdb)
    log=RunWI(com)
    for line in log:
        print line,
    print '*********'
    com=' check \n fulchk \n %s \n \n \n Y \n' %newpdb
    log=RunWI(com)
    record=False
    scores={}
    for line in log:
        print line,
        if len(line.strip())==0:
            continue
        elif line.split()[0]=='RMS':
            continue
        if line[:6]=='SCRIPT'and record:
            break
        if record:
            sp=line.split(':')
            scores[sp[0].strip()]=float(sp[1].split()[0])
        if line.strip()=='Structure Z-scores, positive is better than average:':
            record=True
    os.chdir(topdir)
    shutil.rmtree(tdir)
    if scores=={}:
        raise Exception('WHATCHECK failed')
    return scores
                    
    

def get_abnormal_phipsi(pdbfile):
    import string
    com=' getmol '+pdbfile+' \n Y \n \n check \n chichk \n'
    log=RunWI(com)
    x=0
    while 1:
        #print string.strip(log[x])
        if string.strip(log[x])!=0:
            if string.strip(log[x])=='phi-psi combinations.':
                x=x+2
                break
        if x>=len(log)-1:
            return {}
        x=x+1
    residues={}
    while 1:
        if string.strip(log[x])!='SCRIPT>>>>':
            log[x]=string.replace(log[x],'(',' ')
            log[x]=string.replace(log[x],')',' ')
            tmp=string.split(log[x])
            restyp=tmp[1]
            resnum=tmp[2]
            chainid=tmp[3]
            if len(chainid)!=1:
                chainid=''
            uniqueid=chainid+':'+string.zfill(resnum,4)
            #print uniqueid
            if residues.has_key(uniqueid):
                raise 'uniqueid is not unique'
            else:
                residues[uniqueid]=restyp
        else:
            break
        if x>=len(log)-1:
            break
        x=x+1
    return residues

def get_abnormal_rotamers(pdbfile):
    #
    # Runs WHATIF/CHECK/ROTCHK
    # and returns a list of residues that were found to have unusual rotamers
    #
    import string
    com=' getmol '+pdbfile+' \n Y \n \n check \n rotchk \n'
    log=RunWI(com)
    x=0
    while 1:
        #print string.strip(log[x])
        if len(string.strip(log[x]))!=0:
            #print log[x],
            if string.strip(log[x])=='a summary of the results will be presented.':
                x=x+3
                break
        if x>=len(log)-1:
            return {}
        x=x+1
    residues={}
    while 1:
        if string.strip(log[x])[:3]!='The':
            log[x]=string.replace(log[x],'(',' ')
            log[x]=string.replace(log[x],')',' ')
            tmp=string.split(log[x])
            restyp=tmp[1]
            resnum=tmp[2]
            chainid=''
            if len(tmp[3])==1:
                if tmp[3] in string.letters:
                    chainid=tmp[3]
            try:
                hits=string.atoi(tmp[-1])
            except ValueError:
                break
            #
            # Do we have a score?
            #
            expected=6
            if chainid=='':
                expected=expected-1
            #print tmp
            if len(tmp)==expected:
                try:
                    score=string.atof(tmp[-2])
                except ValueError:
                    score=None
            else:
                score=None
            #
            #
            #
            uniqueid=chainid+':'+string.zfill(resnum,4)
            if residues.has_key(uniqueid):
                pass
            else:
                residues[uniqueid]={'type':restyp,'score':score,'hits':hits}
        else:
            break
        if x>=len(log)-1:
            break
        x=x+1
    return residues    

def get_bumps(pdbfile):
    #
    # Count the number of bumps in the pdbfile
    #
    import os,string
    command='getmol '+os.path.split(pdbfile)[1]+' \n \n set \n %bmpchk \n '
    logfile=RunWI(command,None,None,[pdbfile])
    x=0
    while 1:
        if string.strip(logfile[x])=='one direction.':
            x=x+2
            break
        x=x+1
        if x>=len(logfile)-1:
            return None
    bump=0
    while string.strip(logfile[x])!='SCRIPT>>>>':
        bump=bump+1
        x=x+1
        #print logfile[x],
        if string.lower(string.strip(logfile[x])[:len('and so on for a total of')])=='and so on for a total of':
            return string.atoi(string.split(logfile[x])[7])
    return bump

def get_hbondscore(pdbfile):
    #
    # Gets the hbond score found by HB2NET (given by HB2LFR)
    #
    import os,string
    command='getmol '+os.path.split(pdbfile)[1]+' \n \n set \n setwif 339 1 \n %hb2net \n %hb2lfr'
    logfile=RunWI(command,None,None,[pdbfile])
    x=0
    while 1:
        if string.strip(logfile[x][:len('Total current value....................')])=='Total current value....................':
            split=string.split(logfile[x])
            return string.atof(split[3])
        
        x=x+1
        if x>=len(logfile)-1:
            return None
    return None

def get_saltbridges(pdbfile):
    import os, string
    command=' getmol '+os.path.split(pdbfile)[1]+' \n Y \n set \n symtry \n symspg \n sympar 3\n soushl \n %shosbr \n tot \n tot \n y \n 4.5 \n'
    logfile=RunWI(command,None,None,[pdbfile])
    record=None
    saltbr={}
    for line in logfile:
        if record==1 and string.strip(line)=='SCRIPT>>>>':
            break
        if record==1:
            l2=string.replace(line,'(',' ')
            l2=string.replace(l2,')',' ')
            llist=string.split(l2)
            resnum1=string.zfill(llist[3],4)
            resnam1=llist[2]
            resnum2=string.zfill(llist[8],4)
            resnam2=llist[7]
            if not saltbr.has_key(resnum1):
                saltbr[resnum1]=resnum2
            if not saltbr.has_key(resnum2):
                saltbr[resnum2]=resnum1
        if string.strip(line)=='SCRIPT>>>> 4.5':
            record=1
    return saltbr

#
# --------------------------
#

def get_titratable_residues(pdbfile):
    """
    # Call WHAT IF and get the titratable groups
    """
    import os, string
    command=' getmol '+os.path.split(pdbfile)[1]+' \n Y \n set \n %lsteps '
    logfile=RunWI(command,None,None,[pdbfile])
    record=None
    titratable_groups=[]
    commandok=None
    terminal=None
    for line in logfile:
        #
        # Do we have a 'Bridged' statement for a Cys?
        #
        if line.find('Bridged')!=-1:
            line=line.split('Bridged')[0]
        if string.split(line):
            if string.split(line)[0]=='Number':
                record=None
            if record:
                if string.strip(line)[-4:]=='TERM':
                    terminal=string.strip(line)
                else:
                    if terminal:
                        titratable_groups.append(string.strip(line)+' '+terminal)
                        terminal=None
                    else:
                        titratable_groups.append(string.strip(line))
                commandok=1
            if string.strip(line)=='SCRIPT>>>> %lsteps':
                record=1
    if not commandok:
        raise "WHAT IF failed. Could not get the titratable groups"
    return titratable_groups

#
#  ----------------------
#

def corall(pdbfile,renumb=None):
    """
    # Correct all
    """
    newpdbfile=os.path.split(pdbfile)[1]+'.new'
    command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n 1 \n set \n'
    if renumb:
        command=command+' %renumb \n tot \n 1 \n '
    command=command+' %corall \n N\n %makmol \n '+os.path.split(pdbfile)[1]+'\n'+newpdbfile+' \n \n tot 0 \n'
    logfile,files=RunWI(command,[newpdbfile],None,[pdbfile])
    return logfile,files

def corall_delwat(pdbfile,renumb=None):
    #
    # Correct all and delete water molecules
    #
    newpdbfile=os.path.split(pdbfile)[1]+'.new'
    command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n 1 \n set \n'
    if renumb:
        command=command+' %renumb \n tot \n 1 \n '
    command=command+'  %DELWAT \n%corall \n N\n %makmol \n '+os.path.split(pdbfile)[1]+'\n'+newpdbfile+' \n \n tot 0 \n'
    logfile,files=RunWI(command,[newpdbfile],None,[pdbfile])
    return logfile,files

def dellig_delhyd_corall(pdbfile=None,renumb=None,readall=None,setcha=None):
    #
    # Correct all and delete the ligands
    #
    newpdbfile=os.path.split(pdbfile)[1]+'.new'
    command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n 1 \n set \n'
    if readall:
        command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n 100000 \n set \n'
    if renumb:
        command=command+' %renumb \n tot \n 1 \n '
    if setcha:
        command=command+' %setcha \n tot 0 \n A \n' 
    command=command+' %DELLIG \n %DELHYD tot 0 \n %corall \n N\n %makmol \n '+os.path.split(pdbfile)[1]+'\n'+newpdbfile+' \n \n tot 0 \n'
    logfile,files=RunWI(command,[newpdbfile],None,[pdbfile])
    return logfile,files

def delwat_dellig_corall(pdbfile=None,renumb=None,readall=None,setcha=None):
    #
    # Correct all and delete the ligands
    #
    newpdbfile=os.path.split(pdbfile)[1]+'.new'
    command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n 1 \n set \n'
    if readall:
        command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n 100000 \n set \n'
    if renumb:
        command=command+' %renumb \n tot \n 1 \n '
    if setcha:
        command=command+' %setcha \n tot 0 \n A \n' 
    command=command+' %DELWAT \n %DELLIG \n %corall \n N\n %makmol \n '+os.path.split(pdbfile)[1]+'\n'+newpdbfile+' \n \n tot 0 \n'
    logfile,files=RunWI(command,[newpdbfile],None,[pdbfile])
    return logfile,files

def delwat_dellig_delhyd_corall(pdbfile,renumb=None,readall=None):
    newpdbfile=os.path.split(pdbfile)[1]+'.new'
    command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n 1 \n set \n'
    if readall:
        command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n 100000 \n set \n'
    if renumb:
        command=command+' %renumb \n tot \n 1 \n '
    command=command+' %DELWAT \n %DELLIG \n %DELHYD tot 0 \n %corall \n N\n %makmol \n '+os.path.split(pdbfile)[1]+'\n'+newpdbfile+' \n \n tot 0 \n'
    logfile,files=RunWI(command,[newpdbfile],None,[pdbfile])
    return logfile,files

def corall_addpolhyd_dellig_delwat(pdbfile,renumb=None):
    #
    # Performs a corall and adds hydrogens, and saves a pdb file with only
    # polar hydrogens. This one deletes waters and ions afterwards.. 
    # Stewarts little routine
    #
    newpdbfile='clean.pdb'
    command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n set \n'
    if renumb:
        command=command+' %renumb \n tot \n 1 \n '
    command=command+' %delhyd tot 0 \n %corall \n N\n setwif 339 1 \n %addhyd tot 0 \n tot 0\n setwif 105 1 \n %delwat \n %dellig \n %makmol \n '+os.path.split(pdbfile)[1]+'\n'+newpdbfile+' \n \n tot 0 \n'
    logfile,files=RunWI(command,[newpdbfile],None,[pdbfile,'HTOP'])
    return logfile,files

#
# ---
#

def corall_addpolhyd(pdbfile,renumb=None):
    #
    # Performs a corall and adds hydrogens, and saves a pdb file with only
    # polar hydrogens. 
    # Stewarts little routine
    #
    newpdbfile='clean.pdb'
    command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n set \n'
    if renumb:
        command=command+' %renumb \n tot \n 1 \n '
    command=command+' %delhyd tot 0 \n %corall \n N\n setwif 339 1 \n %addhyd tot 0 \n tot 0\n setwif 105 1 \n %makmol \n '+os.path.split(pdbfile)[1]+'\n'+newpdbfile+' \n \n tot 0 \n'
    logfile,files=RunWI(command,[newpdbfile],None,[pdbfile,'HTOP'])
    return logfile,files

#
# ----
#

def delwat_dellig_corall_addhyd(pdbfile,renumb=None):
    #
    # Deletes all waters and ligand. Performs a corall and adds hydrogens
    # Written for FOLD-X
    #
    newpdbfile=os.path.split(pdbfile)[1]+'.new'
    command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n set \n'
    if renumb:
        command=command+' %renumb \n tot \n 1 \n '
    command=command+' %delhyd tot 0 \n %DELWAT \n %DELLIG \n %corall \n N\n setwif 339 1 \n %addhyd tot 0 \n tot 0\n %makmol \n '+os.path.split(pdbfile)[1]+'\n'+newpdbfile+' \n \n tot 0 \n'
    logfile,files=RunWI(command,[newpdbfile],None,[pdbfile,'HTOP'])
    return logfile,files

#
# ----
#

def delwat_corall_addhyd_setscg(pdbfile,defcrg=0.0,defrad=0.0,renumb=None,newcid='A',crgfile=None,radfile=None):
    #
    # Deletes all waters. Performs a corall and adds hydrogens
    # Written for APBS PDB scale electrostatics
    # UCSD 2001
    #
    newpdbfile=os.path.split(pdbfile)[1]+'.new'
    command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n 1 \n 1 \n'
    if renumb:
        command=command+' %renumb \n tot \n 1 \n %setcha tot 0 \n '+newcid+' \n'
    #
    # 1828 1: Assign default charge and radii for unknown atoms
    # 1829 default charge*100
    # 1830 default radius*100 (in A)
    #
    command=command+'setwif 1828 1\n setwif 1829 '+str(int(100.0*defcrg))+' \n setwif 1830 '+str(int(100.0*defrad))+' \n '
    #
    # 1883 3: Just write a PQR file
    #
    command=command+' %delhyd tot 0 \n %DELWAT \n %corall \n N\n setwif 339 1 \n %addhyd tot 0 \n tot 0\n %setscg \n setwif 1883 3 \n %rundel \n tot 0 \n'
    if crgfile and radfile:
        logfile,files=RunWI(command,['pqr.pdb'],None,[pdbfile,'HTOP',crgfile,radfile])
    else:
        logfile,files=RunWI(command,['pqr.pdb'],None,[pdbfile,'HTOP'])
    return logfile,files

#
# ----

def dellig_addallhyd(pdbfile,renumb=1):
    """Remove all ligands and water, correct all, and add all hydrogens smartly (to correct flips), delete hyds and then add them back stupidly (incl. Asp and Glu)"""
    newpdbfile=os.path.split(pdbfile)[1]+'.new'
    command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n 1 \n 1 \n dellig \n'
    if renumb:
        command=command+' %renumb \n tot \n 1 \n '
    # Addhyd
    command=command+'modpka \n makpsf \n  y \n n \n 1.0 \n \n'
    command=command+'putpsc \n'
    #
    # Save the new PDB file
    #
    command=command+'%makmol \n \n'+newpdbfile+' \n \n tot 0 \n'
    print 'writing',newpdbfile
    #
    # Run WHAT IF
    #
    logfile,files=RunWI(command,[newpdbfile],None,[pdbfile,'HTOP'])
    return logfile,files
    
#

#
# ---------------------------------------------------
#

def delwat_corall_pastal_addhyd_setscg(pdbfile,defcrg=0.0,defrad=0.0,renumb=None,newcid='A',crgfile=None,radfile=None):
    print
    print 'NB: PASTAL HAS BEEN DISABLED!!!'
    print
    #
    # Deletes all waters. Performs a corall and adds hydrogens
    # Written for APBS PDB scale electrostatics
    # UCSD 2001
    #
    newpdbfile=os.path.split(pdbfile)[1]+'.new'
    command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n 1 \n 1 \n'
    if renumb:
        command=command+' %renumb \n tot \n 1 \n %setcha tot 0 \n '+newcid+' \n'
    #
    # 1828 1: Assign default charge and radii for unknown atoms
    # 1829 default charge*100
    # 1830 default radius*100 (in A)
    #
    command=command+'setwif 1828 1\n setwif 1829 '+str(int(100.0*defcrg))+' \n setwif 1830 '+str(int(100.0*defrad))+' \n '
    #
    # 1883 3: Just write a PQR file
    #
    command=command+' %delhyd tot 0 \n %DELWAT \n %corall \n N\n setwif 339 1 \n  %addhyd tot 0 \n tot 0\n %setscg \n setwif 1883 3 \n %rundel \n tot 0 \n'
    if crgfile and radfile:
        logfile,files=RunWI(command,['pqr.pdb'],None,[pdbfile,'HTOP',crgfile,radfile])
    else:
        logfile,files=RunWI(command,['pqr.pdb'],None,[pdbfile,'HTOP'])
    return logfile,files

def equal(pdbfile1,pdbfile2,superpose='N'):
    #
    # Performs an %EQUAL for the two pdbfiles
    #
    import os, string
    pdb1=os.path.split(pdbfile1)[1]
    pdb2=os.path.split(pdbfile2)[1]
    command='getmol '+pdb1+' \n X \n set1 \n getmol '+pdb2+' \n y \n set2 \n'
    #command=command+'%equal \n m1 \n m2 \n Y \n N \n N \n N \n Y \n 1 \n q \n 2 \n w \n 3 \n e \n '
    command=command+'%equal \n m1 \n m2 \n '+superpose+' \n N \n N \n N \n N \n '
    logfile=RunWI(command,None,None,[pdbfile1,pdbfile2])
    alpha=None
    all=None
    res={}
    for line in logfile:
        split=string.split(line)
        if len(split)>=1:
            if split[0]=='Alpha':
                res['alpha']=string.atof(split[4])
            if split[0]=='ALL':
                res['all']=string.atof(split[4])
                break
            if split[0]=='Side':
                res['sidechain']=string.atof(split[5])
            if split[0]=='Backbone':
                res['backbone']=string.atof(split[4]) 
    return res

def has_cys(pdbfile):
    #
    # This function returns 1 if there are any SS bridges in the protein and
    # None if there aren't any.
    #
    import os
    pdb1=os.path.split(pdbfile)[1]
    command='getmol '+pdb1+' \n Y \n set1 \n %shocys '
    logfile=RunWI(command,None,None,[pdbfile])
    for line in logfile:
        if string.strip(line)=='There are no Cys-Cys bridges':
            # No cys bridges...
            return None
    return 1

def relative_accessibility(pdbfile):
    """
    # This option runs DELLIG, INIACC, SETACC and VACACC
    """
    import os
    pdb=os.path.split(pdbfile)[1]
    command='getmol '+pdb+' \n X \n dellig \nset1 \n access \n setacc \n tot 0\n \n vacacc \n N \n tot \n'
    logfile=RunWI(command,None,None,[pdbfile])
    x=0
    resdict={}
    section_found=False
    while 1:
        sline=string.strip(logfile[x])
        if sline.find('vacacc')!=-1:
            section_found=True
        if not section_found:
            x=x+1
            continue
        if len(string.split(sline))>0:
            if string.split(sline)[0]=='Residue:':
                sline=string.replace(sline,'(',' ')
                sline=string.replace(sline,')',' ')
                split=string.split(sline)
                if split[3]=='OXT' or split=='O2':
                    break
                import pKarun.pKa_utility_functions
                residue=pKarun.pKa_utility_functions.getWI_resid3(sline)
                type=split[2]
                resdict[residue]={}
                while 1:
                    x=x+1
                    sline=string.strip(logfile[x])
                    split=string.split(sline)
                    if split[0][0] in string.letters and split[0][:3]!='Phi' and split[0][:4]!='Atom':
                        #
                        # The line below is a mess, it should be rewritten
                        # in a more sensible way to prevent character overflows
                        # in WHAT IF output to break it
                        #
                        resdict[residue][split[0]]={'rel':string.atof(string.split(sline)[-1]),'vac':string.atof(string.split(sline)[-2]),'prot':string.atof(string.split(sline)[4][:4]),'type':type}
                    if len(split)==3 and split[0][0] not in string.letters:
                        resdict[residue]['sum']={'rel':string.atof(string.split(sline)[2]),'vac':string.atof(string.split(sline)[1]),'prot':string.atof(string.split(sline)[0]),'type':type}
                        break
        x=x+1
        if x==len(logfile):
            break
    return resdict
                    

class structure_search:
    #
    # This is a comment
    #
    def __init__(self,debug=None):
        self.hssp='/data/hssp'
        self.pdb='/vriend/hobby1/nielsen/progwork/activesites/pdb/'
        self.debug=debug
        return

    def get_pocket(self,pdbfile):
        FileError="FileError"
        import os, Protool.structure
        #
        # Does the pdbfile exist? - is there a corresponding hssp file?
        #
        p=Protool.structureIO()
        pdbfile=self.long_pdbfilename(pdbfile)
        if not os.path.isfile(os.path.join(self.pdb,pdbfile)) and not os.path.isfile(os.path.join(os.getcwd(),pdbfile)):
            raise FileError,pdbfile
        hsspfile=self.short_pdbfilename(pdbfile)
        if not os.path.isfile(os.path.join(self.hssp,hsspfile+'.hssp')):
            raise FileError,hsspfile
        #
        # Create the POCKET.PDB file and return it with the output
        #
        pocketfile='POCKET.PDB'
        os.putenv('WIFSIZE','0.1')
        command='%stedev 3 \n getmol '+pdbfile+'\n Y \n \n \n %larcav \n '+hsspfile+' \n \n \n \n \n \n \n '
        output,files=RunWI(command,[pocketfile],self.debug,[os.path.join(self.pdb,pdbfile)])
        if files.has_key(pocketfile):
            return output,files[pocketfile]
        else:
            print 'No pocket file'
            return output,None

    def long_pdbfilename(self,pdbfile):
        if pdbfile[-4:]!='.brk' and pdbfile[-4:]!='.pdb':
            return pdbfile+'.brk'
        return pdbfile

    def short_pdbfilename(self,pdbfile):
        if pdbfile[-4:]=='.brk' or pdbfile[-4:]=='.pdb':
            pdbfile=pdbfile[:-4]
        return pdbfile

#
# =============== End of class: structure_search ====================
#


def get_PDB_length(pdbfile,resnum=None):
    import string
    #
    # Count the number of residues in the PDB file
    # This will fail for proteins with non-consequtive
    # residue number and for proteins with more than one
    # chain identifier. Jens 11/10-00
    #
    counter=0
    inputPDB = open(pdbfile).readlines()
    chainid=0
    protein_length=None
    first_residue=None
    for line in inputPDB:
        line = string.split(line)
        if line[0] == 'ATOM' and counter == 0:
            #
            # Check for chain identifier! - Jens 11/10-00
            #
            if len(line[4])==1 and line[4] in string.letters:
                chainid=1
            first_residue=string.atoi(line[4+chainid])
            new_residue = first_residue
            protein_length = 1
            counter = 1
        elif line[0] == 'ATOM' and counter == 1:
            residue=string.atoi(line[4+chainid])
            if resnum:
                if string.atoi(resnum)==residue:
                    return line[3]
            if residue > new_residue:
                protein_length = protein_length + 1
                new_residue = residue
    if resnum:
        return None
    return protein_length, first_residue


#
# ----------------------------------------------------------------------
#

def Remove_line(mutation,pdbfile,newpdbfile,old,new):
    import re,string
    number_mut = mutation[1:-1]
    # How many characters (space or letter) between the name and number of residues #
    numberofspace = "[\s\w]"*(6-len(number_mut))
    KeepAtoms = ['N','HN','CA','C','O','CB','CG','OG1','CG1','CG2','CD1','CD2','CE1','CE2','CZ']
    WhereToStop = {'ALA':6,'GLY':5,'VAL':10,'PHE':len(KeepAtoms),'SER':8}
    input = open("%s"%(pdbfile),'r').readlines()
    newfile = []
    for line in input:
	split = string.split(line)
	if len(split)>0 :
	    if split[0] == 'ATOM':
		if split[3] == old:
		    ii = 0
		    while 1:
			try:
			    eval(split[3+ii])
			    break
			except:
			    ii = ii + 1
		    if split[3+ii] == number_mut:
			for atom in KeepAtoms[:WhereToStop[new]]:
			    if split[2] == atom:
				match1 = re.compile(r"(%s)(%s)(%s)"%(old,numberofspace,number_mut))
				m = match1.search(line)
				if m:
				    #print m.groups()
				    line = match1.sub(r'%s%s%s'%(new,m.groups()[1],m.groups()[2]),line)
				    if split[2] == 'OG1' and old =='THR':
					match2 = re.compile(" OG1 ")
					line = match2.sub(" OG  ",line)
					newfile.append(line)
				    else:
					newfile.append(line)
		    else:
			newfile.append(line)			 
		else:
		    newfile.append(line)		
	    else:
		newfile.append(line)
	else:
	    newfile.append(line)
    # Because the output should be in realines() format we have to save to a temp file
    files = {}
    files[newpdbfile] = newfile
    return files

#
# ===========================================================================
# 

def mutate(pdbfile,mutations,quality=1,silent=None,qualitycheck=None,wiffile=None,chain=None):
    #
    # Make mutations
    # quality==1: Simple
    # quality==2: Use experimental version
    # quality==3: Use experimental version with debumping
    #
    # If multiple mutations are modelled, then we do the simple ones first
    #
    # If there is more than one 'complex' mutation left after doing everything, then we mutate to Ala first and then to the
    # final residue
    #
    if mutations==[] or mutations==[[]] or not mutations:
        raise "No mutations given: ",mutations
    import string, os, sys
    onetothree={'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}
    threetoone={'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}
    returnfiles={}
    oldquality=quality
    orgpdbfile=pdbfile
    tmpfiles=[]
    new=''
    old=''
    #
    # First loop over the sets of mutations
    #
    for mutation1 in mutations:
        pdbfile=orgpdbfile
        #
        # Sort the mutations in this set, so that we make all trivial mutations first
        #
        simple=['DN','ND','EQ','QE','YF','FY','TV','VT']
        mutations_sorted={}
        import string
        for mutation in mutation1:
            mutations_sorted[string.upper(mutation)]=0
        #
        # First we do all -> Ala and -> Gly and -> Pro
        #
        count=0
        for mutation in mutations_sorted.keys():
            if mutation[-1]=='A' or mutation[-1]=='G' or mutation[-1]=='P':
                count=count+1
                if mutations_sorted[mutation]==0:
                    mutations_sorted[mutation]=count
        #
        # Then we do all the simple mutations
        #
        for mutation in mutations_sorted.keys():
            old=mutation[0]
            number=mutation[1:-1]
            new=mutation[-1]
            newold=old+new
            for sim in simple:
                if newold==sim and mutations_sorted[mutation]==0:
                    count=count+1
                    mutations_sorted[mutation]=count
        #
        # And then the small changes Ala -> Ser, etc.
        #
        smallchange=['AS','AC','AT','FY','YF','LI','IL','DE','ED','QN','NQ','DQ','QD','NE','EN']
        for mutation in mutations_sorted.keys():
            old=mutation[0]
            number=mutation[1:-1]
            new=mutation[-1]
            newold=old+new
            for sm in smallchange:
                if newold==sm and mutations_sorted[mutation]==0:
                    count=count+1
                    mutations_sorted[mutation]=count
        #
        #
        # Anything mutated to a small residue comes before big insertions
        #
        small=['S','T','C','V']
        for mutation in mutations_sorted.keys():
            for sm in small:
                if mutation[-1]==sm and mutations_sorted[mutation]==0:
                    count=count+1
                    mutations_sorted[mutation]=count
        #
        # And then the rest
        #
        for mutation in mutations_sorted.keys():
            if mutations_sorted[mutation]==0:
                count=count+1
                mutations_sorted[mutation]=count
        #
        # Load in mutationsort
        #
        mutationsort=[]
        for x in range(1,count+1):
            for mutation in mutations_sorted.keys():
                if mutations_sorted[mutation]==x:
                    mutationsort.append(mutation)
        if not silent:
            print 'Sorted list of mutations: %s' %str(mutationsort)
        #
        #
        #
        # Loop over all mutations in this set
        #
        for mutation in mutationsort:
            if len(mutation)<3 and string.lower(mutation)!='wt':
                print 'Mutations must be given as X#Y, where X is the original residue,',
                print 'Y is the new residue and ## is the number in the PDB file.<BR>'
                print mutation,' failed.<BR>Aborting.<BR>'
                print 'You entered: "%s"<BR>' %mutation
                return None, None, None
            #
            # Initialise
            #
            mutcom='setwif 1012 0 \n'
            #
            # Start checking the mutation
            #
            if string.lower(string.strip(mutation))=='wt':
                quality=4
                mutation=''
            else:
                quality=oldquality
            if not silent:
                print 'Now mutating: %s <BR>' %mutation
            sys.stdout.flush()
            mutation=string.upper(string.strip(mutation))
            old=mutation[0]
            number=str(int(mutation[1:-1]))
            new=mutation[-1]
            newold=old+new
            if (new=='A' and old!='G') or new=='G': 
                quality=1
            for sim in simple:
                if newold==sim:
                    quality=1
            if not silent:
                print '<BR>Using QUALITY %d<BR>' %quality
            oldtype=get_PDB_length(pdbfile,number)
            #
            # Check that we have this residue in the PDB file
            #
            if chain:
                print 'Org residue ID skipped!'
            else:
                if not oldtype:
                    print 'Residue number not found in PDB file: ',number,'<BR>'
                    return None, None, None
                if onetothree.has_key(oldtype):
                    if onetothree[oldtype]!=old:
                        print 'Error: The original residue type is different from what you entered.<BR>'
                        print 'You entered: %s, but the pdbfile says: %s %s.<BR>' %(mutation,oldtype,str(number))
                        return None, None, None
                else:
                    print 'Error: Unrecognised amino acid in PDB file: %s <BR>' %oldtype
                    return None, None, None
            #
            # Check that the new amino acid type exists
            #
            aas=onetothree.values()
            found=None
            for aa in aas:
                if aa==new:
                    found=1
                    break
            if not found:
                print "I don't know this amino acid type (entered as new amino acid) %s %s" %(new,str(number))
                return None, None, None
            #
            # Prepare the mutation command
            #
            # If the new and old residues are identical, then mutate to ALA first...
            #
            if new==old:
                if not chain:
                    mutcom=mutcom+' mutate O'+str(number)+' \n N \n A  \n'
                else:
                    mutcom=mutcom+' mutate U%d \n O %d \n N \n A  \n' %(chain,number)
            #
            # Construct the WHAT IF command string
            #
            if quality==1:
                if not chain:
                    mutcom=mutcom+' mutate O'+str(number)+' \n N \n '+new+' \n'
                else:
                     mutcom=mutcom+' mutate U%d \n %s \n N \n %s \n' %(chain,number,new)
            elif quality==2 or quality==3:
                if not chain:
                    mutcom=mutcom+' mutate O'+str(number)+' \n Y \n '+new+' \n'
                else:
                    mutcom=mutcom+' mutate U%d \n %s \n Y \n %s \n' %(chain,number,new)
            #
            # Everything OK - make the mutation
            #
            newpdbfile=os.path.split(pdbfile)[1]+mutation+'.pdb'
            command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n 1 \n set \n'+mutcom
            if quality==3:
                #
                # Number of debump steps per angle is now in wifpar 73
                #
                # Number of debump steps per angle is 5
                #
                if chain:
                    print 'NOT yet implemented for more than one chain'
                    stop
                command=command+' %PASTML \n setwif 1890 1\n setwif 73 5 \n %DEBHBO N \n O'+str(number)+' \n N \n '
                command=command+' %DEBHBO N \n O'+str(number)+' \n Y \n '
                if quality==4:
                    #
                    # This is to get the wt
                    #
                    mutation='wt'
            #
            # Fill in the rest of the command
            #
            command=command+'%makmol \n \n'+newpdbfile+' \n \n tot 0 \n'
            #
            # Added by Raph to skip WHAT IF if we have a really trivial mutation
            #
            use_whatif = 1
            if not chain:
                list_trivial = ['YF','IV','TS']
                if (new == 'G'and old != 'P'): 
                    files  = Remove_line(mutation,pdbfile,newpdbfile,threetoone[old],threetoone[new])
                    use_whatif = 0
                if (new == 'A' and old != 'P' and old != 'G'): 
                    files  = Remove_line(mutation,pdbfile,newpdbfile,threetoone[old],threetoone[new])
                    use_whatif = 0
                for trivial in list_trivial:
                    if newold == trivial:
                        files  = Remove_line(mutation,pdbfile,newpdbfile,threetoone[old],threetoone[new])
                        use_whatif = 0
            #
            # Use WHAT IF?
            #
            if quality==1:
                rotamer_quality=9999.9 #simple mutations always get perfect score
                #
            else:
                rotamer_quality=None
            if use_whatif == 1:
                #
                # Set the environment
                #
                if os.environ.has_key('LD_LIBRARY_PATH'):
                    ld=os.environ['LD_LIBRARY_PATH']
                else:
                    ld=''
                os.putenv('LD_LIBRARY_PATH','/usr/lib:/lib:/usr/pub/lib:'+ld)
                logfile,files=RunWI(command=command,returnfiles=[newpdbfile],debug=None,neededfiles=[pdbfile],wiffile=wiffile)
                for line in logfile:
                    line=line.strip()
                    #print line
                    if line.find('Mutation score')!=-1:
                        score=float(line.split()[-1])
                        if rotamer_quality:
                            if rotamer_quality>score:
                                rotamer_quality=score
                        else:
                            rotamer_quality=score
            else:
                rotamer_quality=9999.9 # Simple mutation always get perfect rotamer score
            if not rotamer_quality:
                print '\n\n\nWarning: No rotamer quality score for %s\n\n\n' %(str(mutations))
                rotamer_quality=-5.0
            #
            # See if we have a new file
            #
            if len(files.keys())==0:
                #
                # Something went wrong when modelling
                # Try using the simple modelling if that wasn't what we already did...
                #
                for line in logfile:
                    print line,
                stop
                if quality==1:
                    print 'Modelling failed for this residue'
                else:
                    command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n set \n'
                    if not chain:
                        command=command+'mutate O'+str(number)+' \n N \n '+new+' \n'
                    else:
                        command=command+'mutate U%d \n O %d \n N \n %s \n' %(chain,number,new)
                    logfile,files=RunWI(command=command,returnfiles=[newpdbfile],debug=None,neededfiles=[pdbfile],wiffile=wiffile)
                    if len(files.keys())==0:
                        print 'This mutation cannot be modelled.\nWHAT IF output:\n'
                        for line in logfile:
                            print line,
                        return None, None, None
                    else:
                        print 'Mutation could only be modelled using the simple mutation algorithm:',str(number),new
                    if len(files.keys())>0:
                        if files.keys()[0]!=newpdbfile:
                            print 'WHAT IF did not produce a PDB file.<BR>\n'
                            return None, None, None
            #
            # Save the new pdbfile somewhere to reuse it in case there is a new mutaion to be made
            # on top of this one
            #
            import tempfile
            t=tempfile.mktemp()
            fd=open(t,'w')
            for line in files[newpdbfile]:
                fd.write(line)
            fd.close()
            pdbfile=t
            tmpfiles.append(t)
        #
        # Delete all the tempfiles
        #
        for file in tmpfiles:
            if os.path.isfile(file):
                os.unlink(file)
        #
        # Return the output
        #
        files[newpdbfile]=[files[newpdbfile][0]]+['REMARK This mutation was made using quality: '+str(quality)+'\n']+files[newpdbfile][1:]
        returnfiles[string.join(mutation1,'+')]=files[newpdbfile]
    #
    # All sets of mutations are done
    #
    return returnfiles,quality, rotamer_quality

#
# ----------------------------------------------------------------------
#
    
def RunWI(command,returnfiles=None,debug=None,neededfiles=None,wiffile=None):
    #
    # This function runs WHAT IF in a temporary directory
    # returnfiles, and needed files can be specified as arguments
    #
    import tempfile, os
    #
    # Run WHAT IF and return the logfile
    #
    WIERR="WHAT_IF.error"
    top=os.getcwd()
    tmp=tempfile.mkdtemp()
    os.chdir(tmp)
    if debug:
        print tmp
    #
    # Symlink any needed files
    #
    htop=None
    if neededfiles:
        for file in neededfiles:
            if file=='HTOP':
                htop=1
                pass
            else:
                os.symlink(file,os.path.join(tmp,os.path.split(file)[1]))
    #
    # Run normal WHAT IF
    #
    print 'Using normal WHAT IF'
    x=WHAT_IF.Whatif(wiffile)
    #
    # Do we need the hydrogen topology file
    #
    if htop:
        x.getHTOP(tmp)
    x.RedirectOutput('junk.out')
    try:
        x.Run(command,'WIlogfile')
    except WIERR:
        fd=open('WIlogfile')
        output=fd.readlines()
        fd.close()
        import os
        print os.environ
        for line in output:
            print line
        raise WIERR
    fd=open('WIlogfile')
    output=fd.readlines()
    fd.close()
    if returnfiles:
        files={}
        for file in returnfiles:
            if os.path.isfile(file):
                fd=open(file)
                files[file]=fd.readlines()
                fd.close()
        os.chdir(top)
        if not debug:
            pKa_ostools.rmr(tmp)
        return output,files
    else:
        os.chdir(top)
        if not debug:
            pKa_ostools.rmr(tmp)
        return output


#
# The following is Raphs new mutate code (which is basically a slight modification of the above
# mutate routine.
#
def Remove_line(mutation,pdbfile,newpdbfile,old,new):
    import re,string

    number_mut = mutation[1:-1]
    # How many characters (space or letter) between the name and number of residues #
    numberofspace = "[\s\w]"*(6-len(number_mut))
    KeepAtoms = ['N','HN','CA','C','O','CB','CG','OG1','CG1','CG2','CD1','CD2','CE1','CE2','CZ']
    WhereToStop = {'ALA':6,'GLY':5,'VAL':10,'PHE':len(KeepAtoms),'SER':8}
    
    input = open("%s"%(pdbfile),'r').readlines()
    newfile = []
    for line in input:
	split = string.split(line)
	if len(split)>0 :
	    if split[0] == 'ATOM':
		if split[3] == old:
		    ii = 0
		    while 1:
			try:
			    eval(split[3+ii])
			    break
			except:
			    ii = ii + 1
		    if split[3+ii] == number_mut:
			for atom in KeepAtoms[:WhereToStop[new]]:
			    if split[2] == atom:
				match1 = re.compile(r"(%s)(%s)(%s)"%(old,numberofspace,number_mut))
				m = match1.search(line)
				if m:
				    #print m.groups()
				    line = match1.sub(r'%s%s%s'%(new,m.groups()[1],m.groups()[2]),line)
				    if split[2] == 'OG1' and old =='THR':
					match2 = re.compile(" OG1 ")
					line = match2.sub(" OG  ",line)
					newfile.append(line)
				    else:
					newfile.append(line)
		    else:
			newfile.append(line)			 
		else:
		    newfile.append(line)		
	    else:
		newfile.append(line)
	else:
	    newfile.append(line)
	

    # Because the output should be in realines() format we have to save to a temp file
    files = {}
    files[newpdbfile] = newfile
    return files


#----------------------------------------------

def get_hbonds(pdbfile):
    #
    # Returns a dictionary that for every atom gives its hydrogen bond partner
    # If the atom is not present, then it does not make a hydrogen bond
    #
    command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n set \n'
    command=command+' setwif 339 1\n %hb2net \n %hb2lis \n 52 0 \n all 0 \n'
    logfile=RunWI(command,None,None,[pdbfile])
    lines=[]
    record=None
    import string
    for line in logfile:
        split=string.split(line)
        if len(split)>1:
            if split[-1]=='%hb2lis':
                record=1
        if record:
            lines.append(line)
        
    return lines  

#-----------------------------------------------

def run_concrd(pdbfile,structs):
    #
    # Runs CNCPRE and CNCWED for the pdbfile and makes CONCRD produce structs # of structures
    #
    command='getmol '+os.path.split(pdbfile)[1]+' \n Y \n set \n'
    command=command+' concrd \n params \n numout \n '+str(structs)+' \n end \n cncpre tot \n '
    command=command+' cncwed tot \n Y \n'
    logfile,files=RunWI(command,['disco.xtc','dist.pdb'],None,[pdbfile])
    return logfile,files


if __name__=='__main__':
    print relative_accessibility('/home/nielsen/1BVD.pdb')
