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

import pKaTool
import pKarun

def run_WI_pKa_calc(pdbfile,mutated_res,pKarun_params,do_matrix=None):
    """
    # Calculates the desolvation energy, the background interaction energy and all matrix interactions for
    # a single residue (mutated_res) which is given in the form :0000:RES (e.g. :0035:GLU)
    """
    # Create a temp dir
    #
    import tempfile, os, string
    tempfile.tempdir=os.getcwd()
    dirname=tempfile.mkdtemp(dir='/tmp')
    topdir=os.getcwd()
    os.chdir(dirname)
    if os.environ.has_key('HOSTNAME'):
       print 'Running pKa calculations in ',dirname,'on',os.environ['HOSTNAME']
    import sys
    sys.stdout.flush()
    #
    # Get pdbfile, topology.h and the charge files
    #
    # Are the files present in the pwd?
    #
    import os, shutil
    for filename in ['TOPOLOGY.H','DELRAD.DAT','DELCRG.DAT']:
        full_name=os.path.join(topdir,filename)
        destination=os.path.join(dirname,filename)
        if os.path.isfile(full_name):
            shutil.copyfile(full_name,destination)
    #
    # ----
    # Copy the pdb file 
    #
    shutil.copyfile(pdbfile,os.path.join(os.getcwd(),os.path.split(pdbfile)[1]))
    #
    # Find out which group we want to calculate for
    #
    group=None
    import pKarun
    Y=pKarun.pKa_general.pKatools()
    group_list=Y.get_titratable_residues(pdbfile)
    group=[]
    for residue,groupnum in group_list:
        if mutated_res==residue:
            group.append(str(groupnum))
    #
    if len(group)==1:
        #
        # For HIS it seems that we get erratic results for the background interaction energy if the calcs
        # are for only one residue. We therefore add the two previous groups to give more Hbond opt sampling
        #
        if mutated_res.find('HIS')!=-1:
            onlygroup=int(group[0])
            if onlygroup>1:
                group=['%d' %(onlygroup-1)]+group
            if onlygroup>2:
                group=['%d' %(onlygroup-2)]+group
    if group==[]:
        if mutated_res.split(':')[-1]=='CYS':
            # Probably a new SS bridge formed..
            return 'SS-bridge','SS-bridge'
        raise Exception('Could not find WHAT IF group number for mutated residue: %s' %mutated_res)
    #
    # Add this to the parameters
    #            
    pKarun_params['subset']=1
    pKarun_params['groups']=string.join(group,',')
    #
    # Link the pdb file
    #
    pdbfile=os.path.join(topdir,pdbfile)
    #
    # Start the pKa calculation
    #
    X=None
    import pKarun, os
    X=pKarun.pKarun_base.pKarun(os.getcwd(),pdbfile,pKarun_params)
    IO=pKaTool.pKaIO.pKaIO()
    X.clean()
    #
    # Are we do to do a matric calculation or a desolv/backgr calculation?
    #
    if not do_matrix:
        try:
            X.desolv(1)
            desolv=IO.read_desolv(X.desolvfile)
        except:
            desolv=None
        try:
            X.backgr(1)
            backgr=IO.read_backgr(X.backgrfile)
        except:
            backgr=None
    else:
        atomdata=X.matrix(1,style='pKD')
        matrix=IO.read_matrix(X.matrixfile)
        return matrix[mutated_res],atomdata
    #
    # If we cannot calculate one, then warn - this should be changed, and modelling optimised!
    #
    print 'D/B',desolv,backgr
    if desolv is None or backgr is None:
        print '\n\n\nWARNING: Modelling in pdbfile is imperfect - could not estimate intrinsic pKa or interactions\n\n\n'
        return 0.0,0.0
    else:
        #
        # Check that pKa_IO got NTERM and CTERM right
        #
        if not desolv.has_key(mutated_res):
            import pKD_tools
            group_type=pKD_tools.get_titgroup_type_from_titgroup(mutated_res)
            restype=pKD_tools.get_restype_from_titgroup(mutated_res)
            resid=pKD_tools.get_resid_from_res(mutated_res)
            fixed=False
            if group_type=='NTERM':
                newgroup='%s:%s' %(resid,'CTERM')
                if desolv.has_key(newgroup):
                    desolv[mutated_res]=desolv[newgroup]
                    del desolv[newgroup]
                    backgr[mutated_res]=backgr[newgroup]
                    del backgr[newgroup]
                    fixed=True
            elif group_type=='CTERM':
                newgroup='%s:%s' %(resid,'NTERM')
                if desolv.has_key(newgroup):
                    desolv[mutated_res]=desolv[newgroup]
                    del desolv[newgroup]
                    backgr[mutated_res]=backgr[newgroup]
                    del backgr[newgroup]
                    fixed=True
            else:
                pass
            #
            # Did we fix it?
            #
            if not fixed:
                print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print 'Something went wrong in the desolv and backgr calcs - cannot find residue in desolv file'
                print                                                                     
                desolv=None
                backgr=None
                raise Exception('This must not happen')
    #
    # Get the new pKa file
    #
    os.chdir(topdir)
    import pKarun.pKa_ostools as pKa_ostools
    pKa_ostools.rmr(dirname)
    print 'Returning',desolv,backgr
    return desolv, backgr

#
# ----
#

def run_APBS_pKa_calc(pdbfile,mutated_res):
    """
    # Calculates the desolvation energy and the background interaction energy for
    # a single residue (mutated_res) which is given in the form :0000:RES (e.g. :0035:GLU)
    #
    # Create a temp dir
    """
    import tempfile, os, string
    topdir=os.getcwd()
    #
    # Start the pKa calculation
    #
    import pdb2pka_interface
    P2P=pdb2pka_interface.P2pI(topdir,pdbfile,{})
    desolv,backgr=P2P.pdb2pka_desolv_backgr(mutated_res)
    return desolv, backgr

#
# ------------
#
#
# Pre-import the mutate module of Protool so we don't have to
# re-read the rotamer library every time
#
import Protool.mutate
MUT=Protool.mutate.Mutate()

def make_mutation(pdbfile,mutation,topdir=None):
    """Model a mutation"""

    #
    # Adjust name of pdb file
    #
    import os
    pdbfile=os.path.join(os.getcwd(),pdbfile)
    if not os.path.isfile(pdbfile):
        print 'File not found: %s' %pdbfile
        raise Exception()
    #
    # Get the root dir from the pdb file name if it is not specified
    #
    if not topdir:
        split=os.path.split(pdbfile)
        rootdir=split[0]
    else:
        rootdir=topdir
    #
    # Make sure we have a dir created
    #
    import os
    dirname=os.path.join(rootdir,'%s.pdbs' %(os.path.split(pdbfile)[1]))
    if not os.path.exists(dirname):
        if os.path.lexists(dirname):
            os.unlink(dirname)
        os.mkdir(dirname)
    #
    # Is this a real mutation?
    #
    if not mutation:
        print 'Make_mutation got this mutation',mutation
    #
    # Have we produced this PDB file already?
    #
    newfilename=os.path.join(dirname,mutation+'.pdb')
    rotamer_quality_file=os.path.join(dirname,mutation+'_rotq')
    if os.path.isfile(newfilename) and os.path.isfile(rotamer_quality_file):
        #
        # Yes, just pass that filename and the rotamer quality
        #
        fd=open(rotamer_quality_file)
        rotamer_quality=float(fd.readline())
        fd.close()
        return newfilename,rotamer_quality
    else:
        #
        # No :-(, make the mutation
        #
        import Protool, string
        X=Protool.structureIO()
        X.readpdb(pdbfile)
        results={}
        #
        # find the real mutation
        #
        newresID=mutation
        import pKD_tools
        newres=pKD_tools.get_newrestyp_from_mut(mutation)
        oldres=pKD_tools.get_oldrestyp_from_mut(mutation)
        #
        # If this is not a real mutation, then just return the wt pdbfile
        #
        if oldres==newres:
            return pdbfile, 99.9
        #
        # Import Protool and read the pdb file
        #
        import Protool
        P=Protool.structureIO()
        pdbfile=os.path.join(os.getcwd(),pdbfile)
        P.readpdb(pdbfile)
        #
        # Initialise the mutate routines
        #
        MUT.new_PDB(P)
        import pKD_tools
        resid=pKD_tools.get_resid_from_mut(mutation)
        bump_score=MUT.mutate(resid,newres,orgtype=oldres)
        if bump_score is None:
            bump_score=9999999 # set to very high for impossible mutations
        #
        # Write the new file
        #
        P.writepdb(newfilename,style='whatif')
        #
        # Write the rotamer quality file
        #
        fd=open(rotamer_quality_file,'w')
        fd.write('%5.3f\n' %bump_score)
        fd.close()
        print resid,oldres,'->',newres,'bump score',bump_score
    return newfilename,bump_score


#
# -----------------
#

class pKa_info_base:
    #
    # The base class for all classes keeping track of information
    #
    
    def initialise(self):
        import sys
        self.lockfile=self.datafile+'_lock'
        self.data_added=0
        #
        # Get the data array
        #
        import os
        if os.path.isfile(self.datafile) and self.save_file:
            self.lock()
            fd=open(self.datafile)
            import cPickle
            self.data=cPickle.load(fd)
            fd.close()
            self.unlock()
        else:
            self.data={}
        return

    #
    # --------------
    #
    
    def lock(self):
        """
        # Obtain lock on the datafile
        # We wait max. 50 seconds before lock is released!
        """
        return
        import os, time
        locked=None
        count=0
        while not locked or count>50:
            if not os.path.isfile(self.lockfile):
                fd=open(self.lockfile,'w')
                fd.write('locked\n')
                fd.close()
                locked=1
                return
            time.sleep(1)
            print 'Waiting for lock on %s. Count: %d' %(self.lockfile,count)
            count=count+1
            if count>50:
                #self.unlock()
                raise Exception('Could not lock file')
        if not locked:
            print 'I could not obtain lock on: %s' %self.lockfile
            raise Exception()
        return
    #
    # --------------
    #
    
    def unlock(self):
        """
        # Release lock
        """
        return
        import os
        if os.path.isfile(self.lockfile):
            os.unlink(self.lockfile)
        return


    #
    # --------------
    #

    def check_status(self):
        """If we have added more than five datapoints, then store the data on disk"""
        if self.data_added>=5:
            self.update_data()
            self.data_added=0
        return

    #
    # --------------
    #
    
    def __del__(self):
        """Destructor - make sure we save all"""
        self.update_data()
        return

    #
    # --------------
    #

    def update_data(self):
        """Add new data to the file on disk"""
        if self.data_added>0 and self.save_file:
            self.lock()
            import cPickle, os
            if os.path.isfile(self.datafile):
                fd=open(self.datafile)
                old_data=cPickle.load(fd)
                fd.close()
            else:
                old_data={}
            #
            # Put any new data into our dict
            #
            for key1 in old_data.keys():
                if not self.data.has_key(key1):
                    #
                    # copy entire sub-entry
                    #
                    import copy
                    self.data[key1]=copy.deepcopy(old_data[key1])
                else:
                    #
                    # check if we have all subentries
                    #
                    import types, copy
                    if type(old_data[key1]) is types.DictType:
                        for key2 in old_data[key1].keys():
                            if not self.data[key1].has_key(key2):
                                self.data[key1][key2]=copy.deepcopy(old_data[key1][key2])
            #
            # Done, now save the complete dictionary
            #
            tempfile=self.datafile+'_write'
            fd=open(tempfile,'w')
            cPickle.dump(self.data,fd)
            fd.close()
            import shutil
            shutil.move(tempfile,self.datafile)
            self.unlock()
        #
        # Done
        #
        return


    #
    # -----------------
    #

class pKa_info(pKa_info_base):

    def __init__(self,pdbfile,parent,PBEsolver='DelPhi',save_file=True):
        """Instantiate pKa_info"""
        self.pdbfile=pdbfile
        self.PBEsolver=PBEsolver
        self.parent=parent
        self.save_file=save_file
        import os
        self.topdir=os.path.split(pdbfile)[0]
        self.datafile=os.path.join(os.getcwd(),self.pdbfile+'.intpka_data')
        self.initialise()
        self.data_added=0
        return

    #
    # -------------------------
    #

    def get_desolv_and_backgr_mut(self,mutation,residue=None,measuring_mutated_residue=1):
        """
        # Construct the mutation and get all energies
        """
        import string
        if self.parent.mutfile_names.has_key(mutation):
            mutant_pdbfile=self.parent.mutfile_names[mutation]
            print 'Using this mutant PDB file:',mutant_pdbfile
        else:
            mutant_pdbfile,rotamer_quality=make_mutation(self.pdbfile,mutation,self.topdir)
        if not mutation:
            return None
        if not residue:
            import pKD_tools
            residue='%s:%s' %(pKD_tools.get_resid_from_mut(mutation),pKD_tools.get_newrestyp_from_mut(mutation))
        return self.get_desolv_and_backgr(mutant_pdbfile,residue,measuring_mutated_residue)

    #
    # -------------------------
    #

    def get_desolv_and_backgr(self,pdbfilename,residue,measuring_mutated_residue):
        """
        # Get the desolvation energy and the background interaction
        # energy for the residue
        """
        # Do we have to calculate the results?
        #
        #
        # If we are not measuring the mutated residue, then we have to add the context to the
        # residue name
        #
        if not measuring_mutated_residue:
            print 'We are not measuring the mutated reside'
            residue_name='%s@%s' %(residue,pdbfilename)
        else:
            residue_name=residue
        if not self.data.has_key(residue_name):
            #
            # Yes
            # Calculate the background and desolvation
            #
            if self.PBEsolver=='DelPhi':
                self.parent.O._print ('Calculating desolv and backgr for %s %s' %(pdbfilename,residue),text_level=1,tab=0)
                #self.parent.close_stdout()
                desolv,backgr=run_WI_pKa_calc(pdbfilename,residue,self.parent.pKarun_params)
                #self.parent.open_stdout()
            elif self.PBEsolver=='APBS':
                desolv,backgr=run_APBS_pKa_calc(pdbfilename,residue)
            else:
                raise Exception,'Unknown PBEsolver: '+self.PBEsolver
            #
            # Store the result for future use
            #
            if desolv and backgr:
                self.data[residue_name]={'desolv':desolv[residue],'backgr':backgr[residue]}
            else:
                self.data[residue_name]={'desolv':0.0,'backgr':0.0,}
            self.data_added=self.data_added+1
            self.check_status()
        #
        # Return the data
        #
        return self.data[residue_name]['desolv'],self.data[residue_name]['backgr']

#
# -----------------
#

class pKa_dist(pKa_info_base):

    def __init__(self,pdbfile,parent,save_file=True):
        import os
        self.pdbfile=pdbfile
        self.save_file=save_file
        self.datafile=os.path.join(os.getcwd(),pdbfile+'.distances')
        self.topdir=os.path.split(pdbfile)[0]
        self.parent=parent
        self.initialise()
        return

    #
    # ----
    #

    def get_min_dist(self,target_res,mutation):
        #
        # Get the minimum distance between the target residue and
        # the mutated residue
        #
        # Do we know this result already?
        #
        if self.data.has_key(target_res):
            if self.data[target_res].has_key(mutation):
                return self.data[target_res][mutation]
        #
        # No, so we have to calculate it
        # The mutated PDB file might have been created already
        #
        mutfile=None
        if hasattr(self.parent,'mutfile_names'):
            if self.parent.mutfile_names.has_key(mutation):
                mutfile=self.parent.mutfile_names[mutation]
        if not mutfile:
            mutfile,score=make_mutation(self.pdbfile,mutation,self.topdir)
        if not mutfile:
            return None
        import Protool
        X=Protool.structureIO_fast()
        X.readpdb(mutfile)
        #
        # We save a minimum distance for the target residue
        #
        min_dist=9999.9
        for target_atom in X.residues[X.resnum(target_res)]:
            if X.is_backbone(target_atom) or X.is_hydrogen(target_atom):
                continue
            #
            # Loop over all atoms in the mutated residue
            #
            import pKD_tools
            new_resnum=pKD_tools.get_resid_from_mut(mutation) #':'+pKD_tools.get_resnum_from_mut(mutation)
            for mut_atom in X.residues[new_resnum]:
                if X.is_backbone(mut_atom) or X.is_hydrogen(mut_atom):
                    continue
                #
                # Get distance
                #
                distance=X.dist(mut_atom,target_atom)
                if distance<min_dist:
                    min_dist=distance
        #
        # Check the distance in the wt pdb file - we might have removed atoms
        #
        X2=Protool.structureIO_fast()
        X2.readpdb(self.pdbfile)
        for target_atom in X2.residues[X2.resnum(target_res)]:
            if X2.is_backbone(target_atom) or X2.is_hydrogen(target_atom):
                continue
            #
            # Loop over all atoms in the wild type residue
            #
            import pKD_tools
            wt_resnum=pKD_tools.get_resid_from_mut(mutation) #':'+pKD_tools.get_resnum_from_mut(mutation)
            for wt_atom in X2.residues[wt_resnum]:
                if X2.is_backbone(wt_atom) or X2.is_hydrogen(wt_atom):
                    continue
                #
                # Get distance
                #
                distance=X2.dist(wt_atom,target_atom)
                if distance<min_dist:
                    min_dist=distance
                    
        #
        # Save the result
        #
        if not self.data.has_key(target_res):
            self.data[target_res]={}
        self.data[target_res][mutation]=min_dist
        self.data_added=self.data_added+1
        self.check_status()
        #
        # Return the distance
        #
        return min_dist


#
# -----------------
#

class pKa_tabulated(pKa_info_base):
    #
    # Handles all calc + storage of dpKa values for individual point mutations
    #
    def __init__(self,pdbfile,algorithm,use_titration_curves,recalc_intpka,save_file=True):
        """Set the filename for storing the tabulated data"""
        self.pdbfile=pdbfile
        self.save_file=save_file
        import os
        self.datafile=os.path.join(os.getcwd(),self.pdbfile+'.'+algorithm+'.tabdata')
        if use_titration_curves:
            self.datafile=self.datafile+'.titcurves'
        if recalc_intpka:
            self.datafile=self.datafile+'_reintpka_%d' %recalc_intpka       
        self.initialise()
        return

    #
    # ------------------
    #

    def get_dpKa(self,mutation,target,Design_instance):
        """
        # Return the dpKa for the target residue due to mutation
        """
        #
        # If we've already calculated the value then we return it immediately
        #
        import types
        if type(mutation)==types.ListType:
            print
            print 'get_dpka called with wrong argument',mutation
            print 'Was not expecting list'
            raise 'program error'
        if not mutation:
            return 0.0
        if self.data.has_key(mutation):
            if self.data[mutation].has_key(target):
                return self.data[mutation][target]
        #
        # Calculate dpKa value
        #
        if Design_instance.O.level<2:
            text='*********\n**************\nCalculating for tabdata\n\n'
            text=text+'Calculating dpKa for %10s due to mutation %10s ....' %(target,mutation)
            import sys
            sys.stdout.flush()
        #
        # this is dirty - should not call a parent function!
        #
        mut_dpKas=Design_instance.calculate_dpKas_explicit(mutations=[mutation],find_reporter_groups=1)
        #
        if not mut_dpKas.has_key(target):
            print 'dpKa calculation failed for target %s with mutation: %s '%(target,mutation)
            print 'I got this from calculate_dpKas_explicit:'
            residues=mut_dpKas.keys()
            residues.sort()
            print residues
            return None
        #
        # Store the delta pKa value
        #
        if not self.data.has_key(mutation):
            self.data[mutation]={}
        #
        #
        #
        if Design_instance.O.level>2:
            print '%5.2f' %mut_dpKas[target]
        #
        # Get the delta pKa values for all targets
        #
        for this_target in mut_dpKas.keys():
            self.data[mutation][this_target]=mut_dpKas[this_target]
        #
        # 
        #
        if Design_instance.O.level>2:
            print 'Mutation: %15s, target: %15s, dpKa :%4.1f' %(mutation,target,self.data[mutation][target])
        self.data_added=self.data_added+1
        self.check_status()
        return self.data[mutation][target]

    
