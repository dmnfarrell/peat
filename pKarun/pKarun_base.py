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


from pKa_utility_functions import *
#=====================================================
class pKamisc:

    def __init__(self):
        return

    def parse_parameterlist(self,arglist,skip2first=True):
        #
        # This routine parses the parameters given to the pKa routines,
        # sets default parameters and provides some compability with
        # the old syntax
        #
        import string
        #
        # Default parameters
        #
        params={'indi':8,'exdi':80,'ion':0.144,'lowph':0.1,'highph':20.0,'phstep':0.1,'opth':2,'pair':1,'bumps':0,'dessoup':1,'colour':1,'dbcrit':30,'dbdiel':16,'mcrit':1000,'mdiel':60,'yang':0,'debug':0,'symrep':0,'pairene':1,'subset':0,'pka_method':'MC','auto':0,'quickhopt':0,'allow_unknown_atoms':0,'unknown_rad':0,'unknown_crg':0.0,'server':1}
        #
        # Set new values for parameters
        #
        if skip2first:
            parmlist=string.split(string.join(arglist[2:]),'-')
        else:
            parmlist=string.join(arglist).split('-')
        print parmlist
        for parm in parmlist:
            if parm!='':
                argu=string.split(parm)[0]
                if len(string.split(parm))==2:
                    value=string.split(parm)[1]
                    params[argu]=value
                else:
                    params[argu]=1
        #
        # Correct all the stupid artefacts
        #
        if params.has_key('super'):
            del params['super']
            params['opth']=2
        elif params.has_key('noopt'):
            del params['noopt']
            params['opth']=0
        return params

#======================================================================================================
class pKarun:
    
    def __init__(self,dir,pdbfile,params,silent=None,ignore_errors=None,PBEsolver='DelPhi'):
        self.PBEsolver=PBEsolver
        #
        # Set and check basic variables
        #
        import os, string, sys
        if silent:
            oldstdout=sys.stdout
            junk=open('.junkout','w')
            sys.stdout=junk
        #
	self.params={}
        # Make a copy
        self.params=params.copy()
        #
        self.topdir=dir
        self.pdbfile=pdbfile
	if not os.path.isfile(self.pdbfile):
	    raise Exception('PDB file not found: %s' %self.pdbfile)
        if not os.path.isabs(self.pdbfile):
            self.pdbfile=os.path.join(self.topdir,self.pdbfile)
            if not os.path.isfile(self.pdbfile):
                raise Exception('PDBfile is a directory')


        #
	# Default values for all parameters
        #
        defaults={'indi':8,'exdi':80,'ion':0.144,
                  'lowph':2.0,'highph':10.0,
                  'phstep':1.0,'opth':2,'pair':1,
                  'bumps':0,'dessoup':1,'colour':1,
                  'dbcrit':30,'dbdiel':16,'mcrit':1000,
                  'mdiel':60,'yang':0,'debug':0,'symrep':0,
                  'pairene':1,'uhbd':0,'allow_unknown_atoms':0,
                  'unknown_crg':0,'unknown_rad':0,
                  'corall':0,'noBF':0,'subset':0,'groups':None,
                  'group_cutoff':1.0,'selection_rounds':2,
                  'pka_method':'MC',
                  'auto':0,
                  'quickhopt':0,
                  'server':1}

        #
        # The WHAT IF flag is the same as the UHBD flag
        #
        if self.params.has_key('whatif'):
            if self.params['whatif']==1:
                self.params['uhbd']=2
            del self.params['whatif']

        #
        # Set default values for all parameters that we didn't get
        #
 	for key in defaults.keys():
	    if not self.params.has_key(key):
		self.params[key]=defaults[key]

	#	
	# Check that we didn't get a wrong parameter"""
        #
	for key in self.params.keys():
	    if not defaults.has_key(key):
                raise Exception('Wrong parameter: %s' %key)
		

        #---------------------------------------------------------------------------------------------------------------------------
        # Print the Setup
        #
        if self.PBEsolver=='DelPhi':
            print
            print 'Setup:'
            print 'Calculating pKas for   : ',os.path.split(self.pdbfile)[1]
            print 'Protein dielectic      : ',self.params['indi']
            print 'Solvent dielectric     : ',self.params['exdi']
            print 'Ionic strength(M)      : ',self.params['ion']
            print 'pH (start,end,step)    : ',self.params['lowph'],',',self.params['highph'],',',self.params['phstep']
            print 'Proton optimisation    : ',self.params['opth']
            print 'Detailed treatment of pairs: ',self.params['pair']
            print 'Cut-off for strongly interacting pairs (kT): ',self.params['pairene']

            print 'Soup for desolvation (0: all charged/all neutral 1: like in interaction energy calculations.):',self.params['dessoup']
            print 'Yang et al. focussing  :',self.params['yang']
            print 'Alt diel for desolv + backgr             :',self.params['colour']
            if self.params['colour']:
                print 'Alternative dielectric constant          :',self.params['dbdiel']
                print 'When residue colour is above             :',self.params['dbcrit']
            print 'Alternative matrix diel                  :',self.params['mdiel']
            print 'when residue colour is above             :',self.params['mcrit']
            print 'PBE solver: 0:Delphi 1:UHBD 2:WHAT IF    :',self.params['uhbd']
            print 'Method for calculating titration curves  : %s' %(self.params['pka_method'])
            print '---------------------------------------------------------------------------'
            print 'Calculating only for certain groups      : ',self.params['subset']
            print 'Groups selected                          : ',self.params['groups']
            print 'Cutoff for including more groups         : %s' %self.params['group_cutoff']
            print 'Number of selection rounds               : %d' %self.params['selection_rounds']
	print '============================================================='
        #---------------------------------------------------------------------------------------
        # Should we remodel all residues that make a symmetry contact ?
        #
        if self.params['symrep']==1:
            print
            print 'Remodeling all side residues that make symmetry contacts'
            import prepare
            SR=prepare.struct()
            ok=SR.doall(os.path.split(self.pdbfile)[1])
            if ok:
                self.pdbfile=self.pdbfile+'.SYMREP'
            else:
                print 'Using original PDB file'
            print
            print '======================================================'

        #----------------------------------------------------------------------------------------
        # Set all filenames, paths and variables
        #
        self.startcom='getmol '+os.path.split(self.pdbfile)[1]+' xxx \n'
        self.desolvfile=os.path.split(self.pdbfile)[1]+'.DESOLV.DAT'
        self.backgrfile=os.path.split(self.pdbfile)[1]+'.BACKGR.DAT'
        self.matrixfile=os.path.split(self.pdbfile)[1]+'.MATRIX.DAT'
        self.pkafile=os.path.split(self.pdbfile)[1]+'.PKA.DAT'
        self.titcur=os.path.split(self.pdbfile)[1]+'.TITCURV.DAT'
        self.tasklist=[]
        self.flags=''


        #----------------------------------------------------------------------------------------
        # Perform a little preprocessing of some of the flags
        #
        for flag in self.params.keys():
            if flag=='indi' or flag=='mdiel' or flag=='dbdiel' or flag=='exdi':
                self.params[flag]=string.atof(str(self.params[flag]))*10.0
            elif flag=='group_cutoff' or flag=='unknown_rad' or flag=='unkown_crg':
                self.params[flag]=string.atof(str(self.params[flag]))*100.0
            elif flag=='ion' or flag=='pairene':
                self.params[flag]=string.atof(str(self.params[flag]))*1000.0
        #
        # Define the WHAT IF parameter numbers for all the flags that we set.
        #
        parmkey={'indi':1810,'exdi':1811,'yang':1864,'colour':1865,'dbdiel':1867,'mdiel':1868,'dbcrit':1877,'mcrit':1878,'bumps':1858,'pair':1857,'ion':1806,'opth':1855,'dessoup':1863,'debug':1012,'pairene':1856,'uhbd':1883,'group_cutoff':1827,'selection_rounds':1895,'quickhopt':1901,'unknown_rad':1830,'unknown_crg':1829,'allow_unknown_atoms':1828,'server':710}
        #
        # The flags that are processed here should be ignored when we set flags in WHAT IF
        #
        ignore={'lowph':0,'highph':0,'phstep':0,'symrep':0,'corall':0,'noBF':0,'subset':0,'groups':0,'pka_method':0,'auto':0}
        for flag in self.params.keys():
            if ignore.has_key(flag):
                pass
            elif parmkey.has_key(flag):
                self.flags=self.flags+' setwif '+str(parmkey[flag])+' '+str(int(string.atof(str(self.params[flag]))))+' \n '
            else:
                raise Exception('Invalid flag: %s' %flag)

            
        #------------------------------------------------------------------------------
        # Set the flags that never change
        #
        """Resolution in internal map representation"""
        self.flags=self.flags+'setwif 1816 30000 \n '
        """Write everything in soup"""
        self.flags=self.flags+'setwif 105 0 \n '
        """Suppress jokes"""
        self.flags=self.flags+'setwif 1194 1 \n '
        """Make sure that we use 3 letter code for amino acids"""
        self.flags=self.flags+'setwif 1 1\n '


        #-----------------------------------------------------------------------------
        # If we are on Linux then set flag 1846
        #
        if os.environ.has_key('OSTYPE'):
            if string.lower(os.environ['OSTYPE'])=='linux':
                if self.params['uhbd']==1:
                    print 'Using UHBD on Linux with focussing'
                elif self.params['uhbd']==2:
                    print 'Using the WHAT IF PBE on linux'
                else:
                    #self.flags=self.flags+'setwif 1846 1 \n '
                    #print 'Using crappy Linux version of DelPhi'
                    print 'Using new Linux version of DelPhi'
            else:
                print string.lower(os.environ['OSTYPE'])

        #------------------------------------------------------------------------------
        # Should we correct all residues
        #
        if self.params.has_key('corall'):
            if self.params['corall']==1:
                self.startcom=self.startcom+' %corall \n n \n '
        #------------------------------------------------------------------------------
        # Colour all residues according to these rules:
        # 1. Colour according to B-factor
        # 2. If the residue makes crystal contacts, then give it colour 160
        #
        if self.params.has_key('noBF'):
            if self.params['noBF']==1:
                self.startcom=self.startcom+' colour \n colzns\n tot 0 \n 1 \n  colsmc\n tot 0 \n 160 \n %lsteps \n'
            else:
                self.startcom=self.startcom+' colour \n colzns\n tot 0 \n 1 \n colbft\n tot 0 \n colsmc\n tot 0 \n 160 \n %lsteps \n' 
        else:
            self.startcom=self.startcom+' colour \n colzns\n tot 0 \n 1 \n colbft\n tot 0 \n colsmc\n tot 0 \n 160 \n %lsteps \n'
        #------------------------------------------------------------------------------
	self.runstatus(ignore_errors)
        if silent:
            sys.stdout=oldstdout
            junk.close()
            os.unlink('.junkout')
	return

#====================================================================================================
    def runstatus(self,ignore_errors=None):
        import os
        #
        # Find out what has already been done
        #
	self.tasklist=[]
        if not os.path.isfile(self.desolvfile):
            self.tasklist.append('desolv')
        if not os.path.isfile(self.backgrfile):
            self.tasklist.append('backgr')
        if not os.path.isfile(self.matrixfile):
            self.tasklist.append('matrix')
        if not os.path.isfile(self.pkafile) or not os.path.isfile(self.titcur):
            self.tasklist.append('solvepka')
        if len(self.tasklist)==0:
            print
            print '*******************************************************'
            print 'pKa calculations already finished for this PDB file!!!!'
            print '*******************************************************'
        if not os.path.isfile(os.path.join(self.topdir,'TOPOLOGY.H')):
            print 'Please copy TOPOLOGY.H from whatif/dbdata to the pwd'
            print 'pwd is ',self.topdir
            print
            if not ignore_errors and self.PBEsolver=='DelPhi':
                raise Exception('TOPOLOGY.H is not present in the working directory')

        #
        # Write the parameter file
        #
        params_file='run_parameters'
        fd=open(params_file,'w')
        import pickle
        pickle.dump(self.params,fd)
        fd.close()
                                        

        return
    
    def Runall(self):
        #
        # Run everything
        #
	ep=self.reporttime('Program started at:')
	self.all()
	print
	ep2=self.reporttime('Finished at:')
	print 'Time elapsed (mins): %5.2f' %((ep2-ep)/60.0)
	return

    def runseq(self):
	    #
            # Run everything sequentially
            #
            for task in self.tasklist:
                if task=='desolv':
                    """Desolvation"""
                    self.desolv(1)
                elif task=='backgr':
                    """Background interaction"""
                    self.backgr(1)
                elif task=='matrix':
                    """Charge-charge energy matrix"""
                    self.matrix(1)
            #
            # First pass done - now pKa values
            #
            for task in self.tasklist:   
                if task=='solvepka':
                    print 'Running ',task
                    self.solvepka()
                    import os
                    if not os.path.isfile(self.pkafile) or not os.path.isfile(self.titcur):
                        print 'ERROR while running calpka, process:'+str(id)
                        raise Exception('Error while running calpKa, proesss %s' %(str(id)))
            junk=self.reporttime('Time:')
	    return
    
    def clean(self):
        #
        # Remove all resultfiles
        #
	import os
	if os.path.isfile(self.desolvfile):
	    os.unlink(self.desolvfile)
        if os.path.isfile(self.backgrfile):
	    os.unlink(self.backgrfile)
        if os.path.isfile(self.matrixfile):
	    os.unlink(self.matrixfile)
	if os.path.isfile(self.pkafile):
	    os.unlink(self.pkafile)
        if os.path.isfile(self.titcur):
	    os.unlink(self.titcur)
	self.runstatus()
	return

    #
    # ==============================================
    #
    
    def desolv(self,id=0):
        """
        # Calculate desolvation energies
        """
        import os, string
        os.chdir(self.topdir)
        #
        # Did we already do this?
        #
        topresult=os.path.join(self.topdir,self.desolvfile)
        if os.path.isfile(topresult):
            return self.desolvfile
        #
        # Run the command
        #
        self.createdir('desolv',id)
        command=self.startcom+' \n pka \n desolv \n'
        if self.params['groups']:
            #
            # If we calculate only for specific groups then put them in
            # command file
            #
            command=command+' n \n'
            groups=string.split(self.params['groups'],',')
            for group in groups:
                command=command+str(group)+' \n'
            command=command+'0 \n'
        else:
            command=command+'y \n'
        command=self.flags+command+self.desolvfile+' \n '
        import WHAT_IF
        x=WHAT_IF.Whatif()
        x.RedirectOutput('desolv.log')
        x.Run(command,'D.log')
        if not os.path.isfile(self.desolvfile):
            raise Exception('ERROR while running desolv, process: %s' %(str(id)))
        os.link(self.desolvfile,os.path.join(self.topdir,self.desolvfile))
        os.chdir(self.topdir)
	print
	print 'Desolvation energies calculated',
	junk=self.reporttime('Time:')
	print
        return 

    #
    # --------------
    #

    def backgr(self,id=1):
        """
        # Calculate background interaction energies
        """
        import os, string
        os.chdir(self.topdir)
        #
        # Did we calculate background interaction energies already?
        #
        topresult=os.path.join(self.topdir,self.backgrfile)
        if os.path.isfile(topresult):
            return self.backgrfile
        #
        # Go ahead and calculate
        #
        self.createdir('backgr',id)
        command=self.startcom+' \n pka \n backgr \n '
        if self.params['groups']:
            #
            # If we calculate only for specific groups then put them in
            # command file
            #
            command=command+' n \n'
            groups=string.split(self.params['groups'],',')
            for group in groups:
                command=command+str(group)+' \n'
            command=command+'0 \n'
        else:
            command=command+'y \n'
        command=self.flags+command+self.backgrfile+' \n '
        import WHAT_IF
        x=WHAT_IF.Whatif()
        x.RedirectOutput('backgr.log')
        x.Run(command,'B.log')
        if not os.path.isfile(self.backgrfile):
            raise Exception('ERROR when running backgr, process %s' %(str(id)))
        os.link(self.backgrfile,os.path.join(self.topdir,self.backgrfile))
        os.chdir(self.topdir)
	print
	print 'Background energies calculated',
	junk=self.reporttime('Time:')

        if id!=1:
            os._exit(0)
        return 

    #
    # ------------
    #

    def matrix(self,id=1,style=''):
        """
        # Calculate interactions between titratable groups
        """
        import os, string
        os.chdir(self.topdir)
        #
        # Did we already calculate this?
        #
        topresult=os.path.join(self.topdir,self.matrixfile)
        if os.path.isfile(topresult):
            return self.matrixfile
        #
        # Calculate the matrix
        #
        self.createdir('matrix',id)
        if style.lower()=='pkd':
            command=self.startcom+' \n setwif 1902 1 \n pka \n elenet \n'
        else:
            command=self.startcom+' \n pka \n elenet \n'
        if self.params['groups']:
            #
            # If we calculate only for specific groups then put them in
            # command file
            #
            command=command+' n \n'
            groups=string.split(self.params['groups'],',')
            for group in groups:
                command=command+str(group)+' \n'
            command=command+'0 \n'
        else:
            command=command+'y \n'
        command=self.flags+command+self.matrixfile+' \n '
        import WHAT_IF
        x=WHAT_IF.Whatif()
        x.RedirectOutput('matrix.log')
        x.Run(command,'M.log')
        if not os.path.isfile(self.matrixfile):
            os.system('tail -30 matrix.log')
            raise Exception('ERROR while running elenet, process: %s' %(str(id)))
        os.link(self.matrixfile,os.path.join(self.topdir,self.matrixfile))
        #
        # Capture the log file if we are called from pKD
        #
        atomdata=None
        if style.lower()=='pkd':
            fd=open('M.log')
            #
            # Parse the lines
            #
            line=fd.readline()
            while 1:
                line=fd.readline()
                if not line:
                    raise Exception('MATRIX file corrupt')
                if string.find(line,'+++++++++++++')!=-1:
                    break
            #
            # Get potentials at all atoms
            #
            atomdata={}
            while 1:
                exclude={'N':1,'C':1,'H':1,'HN':1,'O':1,'OXT':1}
                line=fd.readline()
                if not line:
                    raise Exception('MATRIX file corrupt!')
                if string.find(line,'^^^^^^^^^^^^^^^^')!=-1:
                    break
                atomid2=getWI_resid4(line)
                atom=string.split(atomid2,':')[-1]
                # We only want the chainID:number
                resid2=string.join(string.split(atomid2,':')[:-2],':')
                energy=string.atof(string.split(line)[-1])
                if not atomdata.has_key(resid2):
                    atomdata[resid2]={}
                # no H's
                if atom[0]=='H' or atom[0]=='1' or atom[0]=='2':
                    continue
                if not exclude.has_key(atom):
                    atomdata[resid2][atom]=abs(energy)
                line=fd.readline()
                if not line:
                    raise Ecxception('MATRIX file corrupt')
            fd.close()
        #
        # End of pKD section
        #
        os.chdir(self.topdir)
	print
	print 'Matrix calculated',
	junk=self.reporttime('Time:')
	print
        if id!=1:
            os._exit(0)
       
        return atomdata

    #
    # ----------------
    #

    def solvepka (self):
        #
        # Calculate titration curves
        #
        import os, string
        os.chdir(self.topdir)
        command='getmol '+os.path.split(self.pdbfile)[1]+' xxx \n pka \n calpka \n'
        if self.params['groups']:
            #
            # If we calculate only for specific groups then put them in
            # command file
            #
            command=command+' n \n'
            groups=string.split(self.params['groups'],',')
            for group in groups:
                command=command+str(group)+' \n'
            command=command+'0 \n'
        else:
            command=command+'y \n'
        #
        # Set the filenames and and pHstart, pHend and pHstep values (the first step is always Tanford-Roxby)
        #
        command=command+'%s \n %s \n %s \n %.2f \n %.2f \n %.2f \n' %(self.desolvfile,self.backgrfile,self.matrixfile,
                                                                      self.params['lowph'],self.params['highph'],self.params['phstep'])
        #
        # Do we want to use the Monte Carlo or the cluster method?
        #
        if self.params['pka_method']=='MC':
            command=command+'Y \n %.2f \n %.2f \n %.2f \n' %(self.params['lowph'],self.params['highph'],self.params['phstep'])
        elif self.params['pka_method']=='TR':
            command=command+'N \n N \n'
        elif self.params['pka_method']=='CLUSTER':
            command=command+'N \n Y \n %.2f \n %.2f \n %.2f \n' %(self.params['lowph'],self.params['highph'],self.params['phstep'])
        #
        # Set the output filenames
        #
        command=command+' %s \n %s \n' %(self.pkafile,self.titcur)
        #
        # Add flags
        #
        command=self.flags+command
        #
        # Run the calcs
        #
        import WHAT_IF
        x=WHAT_IF.Whatif()
        x.RedirectOutput('PKA.log')
        x.Run(command,'P.log')
        return

    #
    # ----
    #

    def autonomy (self):
        """
        # Calculate titration curves for the autonomy analysis
        """
        import os, string
        os.chdir(self.topdir)
        #
        # Right now we need a full pKa calculation to do the autonomy calcs
        #
        command=self.startcom+'\n pka \n autpka \n'
        command=command+'y \n'
        command=command+self.desolvfile+' \n'+self.backgrfile+' \n '+self.matrixfile+' \n '
        command=command+str(self.params['lowph'])+' \n '
        command=command+str(self.params['highph'])+' \n'
        command=command+str(self.params['phstep'])+'\n'
        if self.params['groups']:
            #
            # If we calculate only for specific groups then put them in
            # command file
            #
            command=command+' n \n'
            groups=string.split(self.params['groups'],',')
            for group in groups:
                command=command+str(group)+' \n'
            command=command+'0 \n'
        else:
            command=command+'y \n'
        command=self.flags+command
        import WHAT_IF
        x=WHAT_IF.Whatif()
        x.RedirectOutput('CLUSTER.log')
        x.Run(command,'C.log')
        return        

    #
    # -------------------------------------------
    #       



    #
    # -------------------------------------------
    #
    def grpinf2 (self):
        #
        # Calculate titration curves for the autonomy analysis
        #
        import os, string
        os.chdir(self.topdir)
        #
        # Right now we need a full pKa calculation to do the autonomy calcs
        #
        command=self.startcom+'\n pka \n GRPIN2 \n'
        command=command+'y \n'
        command=command+self.desolvfile+' \n'+self.backgrfile+' \n '+self.matrixfile+' \n '
        command=command+str(self.params['lowph'])+' \n '
        command=command+str(self.params['highph'])+' \n'
        command=command+str(self.params['phstep'])+'\n'
        if self.params['groups']:
            #
            # If we calculate only for specific groups then put them in
            # command file
            #
            command=command+' n \n'
            groups=string.split(self.params['groups'],',')
            for group in groups:
                command=command+str(group)+' \n'
            command=command+'0 \n'
        else:
            command=command+'y \n'
        command=self.flags+command
        import WHAT_IF
        x=WHAT_IF.Whatif()
        x.RedirectOutput('GRPIN2.log')
        x.Run(command,'GRP2.log')
        return        
    #
    # -------------------------------------------
    #
    def grpinf(self):
        #
        # Calculate titration curves for the autonomy analysis
        #
        import os, string
        os.chdir(self.topdir)
        #
        # Right now we need a full pKa calculation to do the autonomy calcs
        #
        command=self.startcom+'\n pka \n GRPINF \n'
        command=command+'y \n'
        command=command+self.desolvfile+' \n'+self.backgrfile+' \n '+self.matrixfile+' \n '
        command=command+str(self.params['lowph'])+' \n '
        command=command+str(self.params['highph'])+' \n'
        command=command+str(self.params['phstep'])+'\n'
        if self.params['groups']:
            #
            # If we calculate only for specific groups then put them in
            # command file
            #
            command=command+' n \n'
            groups=string.split(self.params['groups'],',')
            for group in groups:
                command=command+str(group)+' \n'
            command=command+'0 \n'
        else:
            command=command+'y \n'
        command=self.flags+command
        import WHAT_IF
        x=WHAT_IF.Whatif()
        x.RedirectOutput('GRPINF.log')
        x.Run(command,'GRP.log')
        return        
    #
    # -------------------------------------------
    #
    
    def all(self):
        #
        # Run three parallel jobs
        #
        import os
        pid=os.fork()
        if pid==0:
            id=2
        else:
            id=1
            pid1=pid
        if id==1:
            pid=os.fork()
            if pid==0:
                id=3
            else:
                id=1
                pid2=pid
        for task in self.tasklist:

            if task=='desolv' and id==1:
                """Desolvation"""
                self.desolv(id)
                
            elif task=='backgr' and id==2:
                """Background interaction"""
                self.backgr(id)

                
            elif task=='matrix' and id==3:
                """Charge-charge energy matrix"""
                self.matrix(id)

                
            elif task=='solvepka' and id==1:
                print 'solvepka: Waiting for all processes to finish.'
                st1=os.waitpid(pid1,0)
                st2=os.waitpid(pid2,0)
                print
                print 'Children finished - starting to calculate pKa values'
                print
                """Calculate pKa values"""
                print 'Running ',task,' process:',id
                self.solvepka()
                if not os.path.isfile(self.pkafile) or not os.path.isfile(self.titcur):
                    raise Exception('ERROR while running calpka, process %s' %(str(id)))
                
            junk=self.reporttime('Time:')
        print 'Done'
        if pid==0:
            os._exit(0)
        return

    def createdir(self,task,id):
        #
        # Create a directory
        #
        import os
        print 'Running ',task,' process:',id
        if os.path.isdir(os.path.join(self.topdir,task)):
            os.chdir(task)
            for file in os.listdir('.'):
                os.unlink(file)
            os.chdir(self.topdir)
            os.rmdir(task)
        #
        os.mkdir(task)
        os.chdir(task)
	for file in [self.pdbfile,'../DELRAD.DAT','../DELCRG.DAT','../TOPOLOGY.H','../TOPOLOGY.FIL']:
            if os.path.isfile(file) and not os.path.isfile(os.path.join(os.getcwd(),os.path.split(file)[1])):
                os.symlink(file,os.path.join(os.getcwd(),os.path.split(file)[1]))
        return


    def reporttime(self,text):
        #
        # Print time
        #
        import string, time
        tid=time.localtime(time.time())
        ep=time.time()
        if len(str(tid[4]))==1:
            tt=str(tid[2])+'/'+str(tid[1])+'-'+str(tid[0])+' '+str(tid[3])+':0'+str(tid[4])
        else:
            tt=str(tid[2])+'/'+str(tid[1])+'-'+str(tid[0])+' '+str(tid[3])+':'+str(tid[4])
        print text,tt
        return ep

    def runactsite(self,all_groups):
        #
        # Find the groups that interact strongly with the preselected set
        # of groups
        #
        print 'Finding groups that interact strongly'
        print 'Starting groups:'
        import string
        if not self.params['groups']:
            return None
        groups_tmp=string.split(self.params['groups'],',')
        groups=[]
        for group in groups_tmp:
            groups.append(string.atoi(group))
        #
        for group in groups:
            found=None
            for uniqueid,groupnum in all_groups:
                if groupnum==group:
                    print '%12s (#%4d)' %(uniqueid,groupnum)
                    found=1
                    break
            if not found:
                raise Exception('Invalid group number for this PDB file: %s' %(str(group)))
        import os
        os.chdir(self.topdir)
        #
        # Prepare the command
        #
        command=''
        command=command+self.startcom+' \n pka \n'
        command=command+'actpka \n '
        for group in groups:
            command=command+str(group)+' \n'
        command=command+'0 \n '
        command=self.flags+command
        #
        # Run the command
        #
        print
        print 'Now running the calculation... please wait'
        import WHAT_IF
        x=WHAT_IF.Whatif()
        x.RedirectOutput('groups.log')
        x.Run(command,'G.log') 
        fd=open('G.log')
        lines=fd.readlines()
        fd.close()
        count=0
        groups=0
        print 'Done.\n'
        import string
        #
        # Look for the line where all the groups are listed
        #
        for line in lines:
            split=string.split(line,':')
            if len(split)==2:
                if split[0]=='Number of important groups':
                    groups=string.atoi(split[1])
                    count=count+1
                    break
            count=count+1
        #
        # Get the groups
        #
        group_list=[]
        for x in range(groups):
            group=getWI_resid(lines[count])
            group_list.append(group)
            count=count+1
        #
        # Check that there are no duplicates
        #
        group_list.sort()
        newlist=[]
        for x in range(len(group_list)-1):
            thisgroup=group_list[x]
            found=None
            for testgroup in group_list[x+1:]:
                if testgroup==thisgroup:
                    found=1
            if not found:
                newlist.append(thisgroup)
        newlist.append(group_list[-1])

        #
        # Return
        #
        return newlist
    #
    # ---------------------------------
    #

    def sugelm(self,id=1):
        #
        # Calculate effects of mutations
        #
        import os, string
        os.chdir(self.topdir)
        if not os.path.isdir('sugelm'):
            self.createdir('sugelm',id)
            #
            # We always use the quick Hopt here!
            #
            command=self.startcom+' \n pka \n sugelm \n'
            if self.params['groups']:
                #
                # If we calculate only for specific groups then put them in
                # command file
                #
                command=command+' n \n'
                groups=string.split(self.params['groups'],',')
                for group in groups:
                    command=command+str(group)+' \n'
                command=command+'0 \n'
            else:
                command=command+'y \n'
            command=self.flags+command
            import WHAT_IF
            x=WHAT_IF.Whatif()
            x.RedirectOutput('sugelm.log')
            x.Run(command,'S.log')
            print 'SUGELM finished'
        else:
            os.chdir('sugelm')
        #
        # Parse the output file
        #
        print 'Parsing SUGELM log file'
        data={}
        atomdata={}
        fd=open('S.log')
        line=fd.readline()
        while line:
            if string.find(line,'&&&&&&&&&&')!=-1:
                #
                # Get the number and name of the residue that was mutated
                #
                line=fd.readline()
                resid=getWI_resid(line)
                if not data.has_key(resid):
                    data[resid]={'Rotamer quality':-5.0}
                    atomdata[resid]={'Rotamer quality':-5.0}
                #
                # Get all the interaction energies
                #
                while 1:
                    nterm=None
                    cterm=None
                    line=fd.readline()
                    if string.find(line,'Reading map with')!=-1:
                        break
                    if string.find(line,'NTERM')!=-1:
                        nterm=1
                        line=fd.readline()
                    if string.find(line,'CTERM')!=-1:
                        cterm=1
                        line=fd.readline()
                    resid2=getWI_resid2(line)
                    energy=string.atof(string.split(line)[-1])
                    if nterm:
                        resid2=resid2+':NTERM'
                    if cterm:
                        resid2=resid2+':CTERM'
                    data[resid][resid2]=energy
                #
                # Scan to next section
                #
                while 1:
                    line=fd.readline()
                    if string.find(line,'+++++++++++++')!=-1:
                        break
                #
                # Get potentials at all atoms
                #
                while 1:
                    exclude={'N':1,'C':1,'H':1,'HN':1,'O':1,'OXT':1}
                    line=fd.readline()
                    if string.find(line,'^^^^^^^^^^^^^^^^')!=-1:
                        break
                    atomid2=getWI_resid4(line)
                    atom=string.split(atomid2,':')[-1]
                    # We only want the chainID:number
                    resid2=string.join(string.split(atomid2,':')[:-2],':')
                    energy=string.atof(string.split(line)[-1])
                    if not atomdata[resid].has_key(resid2):
                        atomdata[resid][resid2]={}
                    # no H's
                    if atom[0]=='H' or atom[0]=='1' or atom[0]=='2':
                        continue
                    if not exclude.has_key(atom):
                        atomdata[resid][resid2][atom]=abs(energy)
            line=fd.readline()
        fd.close()
        #
        # Check again for errors
        #
        fd=open('S.log')
        line=fd.readline()
        while line:
            if line.find('Some error from CALCINT')!=-1:
                fd.close()
                fd=open('S.log')
                line=fd.readline()
                while line:
                    line=fd.readline()
                    print line,
                fd.close()
                raise Exception('ERROR when running SUGELM')
            line=fd.readline()
        #
        # Finish
        #
        os.chdir(self.topdir)
	print
	print 'Sugelm finished',
	junk=self.reporttime('Time:')
        print 'Totally done'
        return data,atomdata
