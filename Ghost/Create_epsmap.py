#!/usr/bin/env python
#
# Create a 3D epsmap
#
import numpy
import sys, os
topdir=os.getcwd()
#
# Cluster paths
#
APBSpath='/home/nielsen/bin/APBS_install/bin/apbs'
VGRIDpath='/home/nielsen/bin/APBS/tools/python/vgrid'

if os.environ.has_key('HOST'):
    if os.environ['HOST']=='amylase':
        APBSpath='/home/nielsen/bin/APBS_install/bin/apbs'
        VGRIDpath='/home/nielsen/bin/APBS/tools/python/vgrid'
    elif os.environ['HOST']=='predrag':
        APBSpath='predrag'
        VGRIDpath='predrag'
if os.environ.has_key('OSTYPE'):
    if os.environ['OSTYPE']=='darwin': #Jens' Macbook Pro
        APBSpath='/Users/nielsen/bin/APBS_install/bin/apbs'
        VGRIDpath='/Users/nielsen/bin/APBS/tools/python/vgrid'

sys.path.append(VGRIDpath)

try:
	from vgrid import *
	startVio()
except:
	print 'VGrid routines not available'
	
import sys
from sys import stdout, stderr

#
# Constants
#
ppm_au=1.9447E-18
NSP={'N':977,'H':90,'HA':60}

#
# --------------
#

def length(vector):
    # This function returns the length of vector
    import math
    sum=0.0
    for value in vector:
        sum=sum+math.pow(value,2)
    return math.sqrt(sum)

#
# ------
#
import LM_functions, epsmap2
class iterative_epsmap(LM_functions.LM_functions,epsmap2.epsmap2):

    def __init__(self,eps,options,init_epsmap=False):
        """Change epsmap to get best agreement with experimental Ghosts"""
        #
        # Cluster?
        #
        self.myid=0
        if options.cluster:
            import pypar
            self.myid=pypar.rank()
            self.procs=pypar.size()
            self.node=pypar.get_processor_name()
        #
        # PBE or Coulomb?
        #
        #if options.pbe:
        #    self.method='APBS'
        #else:
        #    self.method='Coulomb'
        self.eps=eps
        self.options=options
        #
        # Read experimental restraints
        #
        
        import os
        self.topdir=topdir
        self.startdir=os.path.join(self.topdir,'start_values')
        if not self.options.make_expdata:
            self.read_experimental_restraints()
        else:
            self.titgroups={}
        #
        # Get the PDB file and clean it up
        #
        import Protool
        self.PI=Protool.structureIO()
        self.PI.readpdb(os.path.join(self.topdir,options.pdb))
        self.PI.RemoveALT()
        self.PI.Remove_All_NonAminoAcids()
        self.PI.setcha(self.PI.residues.keys(),'')
        self.PDB_titgroups=self.PI.get_titratable_groups()
        #
        # Check that the PDB file has all the titgroups
        #
        fail=False
        if self.options.verbose:
            print 'Finding titgroups in PDB file'
        for tg in self.titgroups:
            if not tg in self.PDB_titgroups:
                print 'Could not find %s in PDB file' %tg
                fail=True
        if fail:
            print 'Titgroups in PDB file:'
            print sorted(self.PDB_titgroups)
            #print sorted(self.PI.residues.keys())
            print options.pdb
            raise Exception('Input error')
        #
        # Initialize epsmap
        #
        self.epsmap=None
        if self.options.pbe and self.options.epsmap and not init_epsmap:
            import Epsmap_class
            self.epsmap=Epsmap_class.Epsmap_class(os.path.join(self.startdir,'xdiel.dx'),
                                                    os.path.join(self.startdir,'ydiel.dx'),
                                                    os.path.join(self.startdir,'zdiel.dx'))
        #
        # Done with init
        #
        return
    
    #
    # ------
    #
    
    def read_experimental_restraints(self):
        """Read the experimental Ghosts"""
        import os
        self.expfile=os.path.join(self.topdir,self.options.restraints)
        if self.options.verbose:
            print 'Reading experimental Ghosts from file: :%s' %self.expfile
        if False:
            #
            # This is the old code for reading restraints
            #
            fd=open(self.expfile)
            lines=fd.readlines()
            fd.close()
            #
            # Format of file: titgroup residue N_ghost C_ghost
            #
            self.exp_ghosts={}
            if self.options.verbose:
                print 'Experimental Ghosts'
                print 'Titgroup      Residue       N_ghost N_error   H_ghost  H_Error'
            for line in lines:
                if line[0]=='#':
                    continue
                sp=line.strip().split()
                titgroup=sp[0]
                if not self.exp_ghosts.has_key(titgroup):
                    self.exp_ghosts[titgroup]={}
                residue=sp[1]
                self.exp_ghosts[titgroup][residue]={}
                if sp[2].lower()!='none':
                    self.exp_ghosts[titgroup][residue]['N']=float(sp[2])
                    if sp[3].lower()!='none':
                        self.exp_ghosts[titgroup][residue]['N_error']=float(sp[3])
                    else:
                        self.exp_ghosts[titgroup][residue]['N_error']=0.05
                #
                if sp[4].lower()!='none':
                    self.exp_ghosts[titgroup][residue]['H']=float(sp[4])
                    if sp[5].lower()!='none':
                        self.exp_ghosts[titgroup][residue]['H_error']=float(sp[5])
                    else:
                        self.exp_ghosts[titgroup][residue]['H_error']=0.05
                #
                if self.options.verbose:
                    print '%10s %10s %10s' %(titgroup,residue,self.exp_ghosts[titgroup][sp[1]])
        #
        # Read Predrag's new restraints
        #
        fd=open(self.expfile)
        import pickle
        self.exp_ghosts=pickle.load(fd)
        fd.close()
        #
        # Set all unknown errors to 0.05
        #
        for titgroup in self.exp_ghosts.keys():
            for residue in self.exp_ghosts[titgroup].keys():
                if not self.exp_ghosts.has_key('N_error'):
                    self.exp_ghosts[titgroup][residue]['N_error']=0.05
                if not self.exp_ghosts.has_key('H_error'):
                    self.exp_ghosts[titgroup][residue]['H_error']=0.05
        #
        # Get the names of the titgroups in the experimental file
        #
        self.titgroups=self.exp_ghosts.keys()
        self.titgroups.sort()
        if self.options.verbose:
            print
            print 'Experimental ghosts have been read'
            print 
        return

    #
    # ----
    #

    def slave_wait(self):
        """Wait for messages from the master"""
        while True:
            import pypar
            source=0
            command=pypar.receive(source)
            self.cube_eps=command[0].copy()
            cube=command[1]
            if cube<0:
                print 'Slave %d exiting' %self.myid
                import sys
                sys.exit(0)
            score=self.calc_cube(cube)
            pypar.send(score,0)
        return
    
    #
    # -----
    #

    def calc_cube(self,centercube,master=False):
        """Calculate the score when performing a permutation on this cube"""
        #
        # First init the epsmap
        #
        if not master:
            cubes=self.epsmap.get_cubes(sidelength=self.options.sidelength,eps=80.0)
        for cube in self.cube_eps.keys():
            self.epsmap.set_cubeeps(cube,self.cube_eps[cube])
        #
        # Do the permutation and calculate
        #
        old_diel=self.cube_eps[centercube]
        self.epsmap.set_cubeeps(centercube,old_diel/2.0)
        score=self.get_spans(epsmap=self.epsmap)
        self.epsmap.set_cubeeps(centercube,old_diel)
        return score
        
    #
    # -----
    #
    
    def find_epsmap(self):
        """Find the epsmap"""
        #
        # Make the resultdir
        #
        import os
        #
        # Define the epsmap
        #
        cubes=self.epsmap.get_cubes(sidelength=self.options.sidelength,eps=80.0,clearmaps=False)
        print 'Creating map with %d cubes' %(len(cubes))
        self.cube_eps={}
        self.cube_grad={}
        for cube in cubes:
            self.cube_eps[cube]=80.0
            self.cube_grad[cube]=0.0
            self.epsmap.set_cubeeps(cube,80.0)
        #
        # If starteps is defined, then load the starting epsmap
        #
        if options.starteps=='':
            print 'Using epsmap with 80 everywhere as starting point for iteration'
        else:
            if options.starteps=='PBE':
                diels=[]
                print 'Using normal eps/80 PBE dielectric map as starting EPSMAP'
                count=0
                for cube in cubes:
                    centre=self.epsmap.cubes[cube]['coord']
                    diel=self.epsmap.maps['x'].get_coord_value(centre[0],centre[1],centre[2])
                    if int(diel)!=80:
                        count=count+1
                        diels.append(diel)
                    self.cube_eps[cube]=diel
                print '%4d (%5.2f%%) cubes have a diel different from 80. Their avg diel is %5.2f' %(count,float(count)/float(len(cubes))*100.0,sum(diels)/len(diels))
            else:
                print 'Loading starting epsmap from: %s' %options.starteps
                fd=open(options.starteps)
                import pickle
                self.epsmap=pickle.load(fd)
                fd.close()
                if len(self.epsmap.keys())!=len(cubes):
                    print 'Loaded starting epsmap has wrong size'
                    print len(self.epsmap.keys()),len(cubes)
                    raise Excepion('invalid starting epsmap')
        #
        # Get the starting score
        #
        for cube in cubes:
            self.epsmap.set_cubeeps(cube,self.cube_eps[cube])
        best_score=self.get_score(self.epsmap)
        #
        print 'Original score is: %5.2f' %best_score
        orgscore=best_score
        #
        import os
        converged=False
        #
        # Set starting map
        #
        scores=[]
        diff=1000
        step=0
        import os
        logfile=os.path.join(self.topdir,'mini_%d.log' %self.myid)
        fd=open(logfile,'w')
        fd.write('Starting minimisation. Staring score: %5.3f\n' %orgscore)
        fd.close()
        while diff>0.5:
            #
            # Minimize the system
            #
            cube_scores=[]
            #
            # distribute work to all nodes
            #
            if options.cluster:
                self.cube_grad=self.get_cube_scores()
                for cube in self.cube_grad.keys():
                    self.cube_grad[cube]=self.cube_grad[cube]-best_score
            else:
                #
                # do everything with one prcessor
                #
                for cube in cubes:
                    print 'Getting gradient for cube: %d' %cube
                    old_diel=self.cube_eps[cube]
                    self.epsmap.set_cubeeps(cube,old_diel/2.0)
                    score=self.get_score(self.epsmap)-best_score
                    self.cube_grad[cube]=score
                    self.epsmap.set_cubeeps(cube,old_diel)
            #
            # Change the epsmap
            #
            diff=0.0
            for cube in cubes:
                if abs(self.cube_grad[cube])>0.00001:
                    thisgrad=self.cube_grad[cube]
                    oldeps=self.cube_eps[cube]
                    change=float(oldeps)-oldeps/2.0
                    
                    epschange=thisgrad*float(options.stepsize)
                    epschange=max(-1.0,min(1.0,epschange)) # Make sure we are between 1 and -1
                    epschange=float(epschange*change)
                    #
                    neweps=oldeps+epschange
                    print 'Changing eps from %5.2f to %5.2f with gradient: %5.3f, epschange: %6.3f' %(oldeps,neweps,thisgrad,epschange)
                    diff=diff+abs(epschange)
                    self.cube_eps[cube]=neweps
            #
            # Score this epsmap
            #
            for cube in cubes:
                self.epsmap.set_cubeeps(cube,self.cube_eps[cube])
            score=self.get_score(self.epsmap)
            best_score=score
            #
            # Save this epsmap
            #
            import os
            trajdir=os.path.join(self.topdir,'trajectory')
            if not os.path.isdir(trajdir):
                os.mkdir(trajdir)
            fd=open(os.path.join(trajdir,'map_step_%d.pickle' %step),'w')
            import pickle
            pickle.dump({'sidelength':options.sidelength,'cubeeps':self.cube_eps,'PDBfile':options.pdb},fd)
            fd.close()
            #
            # Save the cube scores
            #
            import os
            cubedir=os.path.join(self.topdir,'cubescores')
            if not os.path.isdir(cubedir):
                os.mkdir(cubedir)
            fd=open(os.path.join(cubedir,'scores_step_%d.pickle' %step),'w')
            import pickle
            pickle.dump(self.cube_grad,fd)
            fd.close()
            #
            step=step+1
            #
            # Append to the info file
            #
            fd=open(logfile,'a')
            fd.write('Step: %d, Score: %5.3f, diff: %5.3f\n' %(step,score,diff))
            fd.flush()
            fd.close()
            print 'Step: %d, Score: %5.3f, diff: %5.3f\n' %(step,score,diff)
        #
        # All done
        #
        print 'Minimisation criterion reached'
        return
        
    #
    # ------
    #
    
    def LM_find_epsmap(self):
        """Find the epsmap"""
        #
        # Make the resultdir
        #
        import os
        logfile=os.path.join(self.topdir,'mini_%d.log' %self.myid)
        fd=open(logfile,'w')
        fd.write('Starting LM minimisation.\n')
        fd.close()
        #
        # Define the epsmap
        #
        cubes=self.epsmap.get_cubes(sidelength=self.options.sidelength,eps=80.0,clearmaps=False)
        print 'Creating map with %d cubes' %(len(cubes))
        self.cube_eps={}
        self.cube_grad={}
        for cube in cubes:
            self.cube_eps[cube]=80.0
            self.cube_grad[cube]=0.0
            self.epsmap.set_cubeeps(cube,80.0)
        #
        # If starteps is defined, then load the starting epsmap
        #
        if options.starteps=='':
            print 'Using epsmap with 80 everywhere as starting point for iteration'
        else:
            if options.starteps=='PBE':
                print 'Using normal eps/80 PBE dielectric map as starting EPSMAP'
                count=0
                diels=[]
                for cube in cubes:
                    centre=self.epsmap.cubes[cube]['coord']
                    diel=self.epsmap.maps['x'].get_coord_value(centre[0],centre[1],centre[2])
                    if int(diel)!=80:
                        count=count+1
                        diels.append(diel)
                    self.cube_eps[cube]=diel
                print '%4d (%5.2f%%) cubes have a diel different from 80. Their average eps is %5.1f' %(count,float(count)/float(len(cubes))*100.0,sum(diels)/len(diels))
            else:
                print 'Loading starting epsmap from: %s' %options.starteps
                fd=open(options.starteps)
                import pickle
                self.epsmap=pickle.load(fd)
                fd.close()
                if len(self.epsmap.keys())!=len(cubes):
                    print 'Loaded starting epsmap has wrong size'
                    print len(self.epsmap.keys()),len(cubes)
                    raise Excepion('invalid starting epsmap')
        #
        # ----
        #
        conv_crit=0.000005
        self.LM_damper = 0.0001
        #
        # Start fitting
        #
        old_diff=0.0000001
        print 'Calculating starting score'
        now_diff=self.get_score()
        step=0
        fd=open(logfile,'a')
        fd.write('Step: %d, diff: %5.3f\n' %(step,now_diff))
        fd.flush()
        fd.close()
        #
        print 'Step    Diff   LM_damper'
        print '%5d  %6.4f  %6.4f' %(0, now_diff, self.LM_damper)
        max_steps=1000
        for step in range(1,max_steps):
            #
            # Change the epsmap
            #
            self.fit_LM_ghost()
            #
            # Calculate new score
            #
            now_diff=self.get_score()
            #
            # Save this epsmap
            #
            import os
            trajdir=os.path.join(self.topdir,'trajectory')
            if not os.path.isdir(trajdir):
                os.mkdir(trajdir)
            fd=open(os.path.join(trajdir,'map_step_%d.pickle' %step),'w')
            import pickle
            pickle.dump({'sidelength':options.sidelength,'cubeeps':self.cube_eps,'PDBfile':options.pdb},fd)
            fd.close()
            #
            # Check convergence
            #
            print '%5d  %6.4f  %6.4f' %(step, now_diff, self.LM_damper)
            if abs(now_diff-old_diff)<conv_crit and False:
                print 'Converged',now_diff
                break
            old_diff=now_diff
            #
            # Save the cube scores
            #
            import os
            cubedir=os.path.join(self.topdir,'cubescores')
            if not os.path.isdir(cubedir):
                os.mkdir(cubedir)
            fd=open(os.path.join(cubedir,'scores_step_%d.pickle' %step),'w')
            import pickle
            pickle.dump(self.cube_grad,fd)
            fd.close()
            #
            #
            # Append to the info file
            #
            fd=open(logfile,'a')
            fd.write('Step: %d, diff: %5.3f\n' %(step,now_diff))
            fd.flush()
            fd.close()
        #
        # All done
        #
        print 'Minimisation criterion reached'
        return
        
    #
    # -----
    #
    
    def cube_scan(self):
        """Calculate the score for all cubes"""
        #
        # Define the epsmap
        #
        cubes=self.epsmap.get_cubes(sidelength=self.options.sidelength,eps=80.0,clearmaps=False)
        print 'Creating map with %d cubes' %(len(cubes))
        self.cube_eps={}
        self.cube_grad={}
        for cube in cubes:
            self.cube_eps[cube]=80.0
            self.cube_grad[cube]=0.0
            self.epsmap.set_cubeeps(cube,80.0)
        #
        # If starteps is defined, then load the starting epsmap
        #
        if options.starteps=='':
            print 'Using epsmap with 80 everywhere as starting point for iteration'
        else:
            if options.starteps=='PBE':
                print 'Using normal 8/80 PBE dielectric map as starting EPSMAP'
                count=0
                diels=[]
                for cube in cubes:
                    centre=self.epsmap.cubes[cube]['coord']
                    diel=self.epsmap.maps['x'].get_coord_value(centre[0],centre[1],centre[2])

                    if int(diel)!=80:
                        count=count+1
                        diels.append(diel)
                    self.cube_eps[cube]=diel
                if len(diels)==0:
                    diels.append(0.0)
                print '%4d (%5.2f%%) cubes have a diel different from 80. Their avg eps is %5.1f' %(count,float(count)/float(len(cubes))*100.0,sum(diels)/len(diels))
                print diels
            else:
                print 'Loading starting epsmap from: %s' %options.starteps
                fd=open(options.starteps)
                import pickle
                self.epsmap=pickle.load(fd)
                fd.close()
                if len(self.epsmap.keys())!=len(cubes):
                    print 'Loaded starting epsmap has wrong size'
                    print len(self.epsmap.keys()),len(cubes)
                    raise Excepion('invalid starting epsmap')
        #
        # ---
        #
        self.cubes=sorted(self.cube_eps.keys())
        self.cube_grad=self.get_cube_scores()
        os.chdir(topdir)
        fd=open('cube_scan_sidelength_%5.1f_eps%5.1f.pickle' %(options.sidelength,options.eps),'w')
        import pickle
        pickle.dump([self.epsmap.cubes,self.cube_grad],fd)
        fd.close()
        #
        # Tell all children to die
        #
        if self.options.cluster:
            import pypar
            print 'Master is killing all children'
            for proc in range(1,self.procs):
                pypar.send([{},-1],proc)
        return
        
    #
    # -----
    #
    
    def get_cube_scores(self):
        """Get scores for all cubes"""
        print
        print 'Starting to calculate scores for all cubes'
        print
        if self.options.cluster:
            import pypar
        else:
            print 'No pypar'
            self.procs=1
            self.ID=0
        print 'Number of procs',self.procs
        proc_status={}
        for proc in range(1,self.procs):
            proc_status[proc]='free'
        #
        # List of cubes
        #
        cubes=sorted(self.cube_eps.keys())
        #
        # Give each free processor a cube and wait for reply
        #
        import sys
        all_done=False
        self.cube_grad={}
        to_get=len(cubes)
        print
        while not all_done:
            txt='Cube scores to get: %4d' %to_get
            print (1+len(txt))*'\b'+txt,
            sys.stdout.flush()
            #
            for proc in range(1,self.procs):
                #print 'proc status',proc_status[proc]
                if proc_status[proc]=='free':
                    for cube in cubes:
                        if not self.cube_grad.has_key(cube) and self.options.cluster:
                            #
                            # No job for this cube yet - send to a node
                            #
                            pypar.send([self.cube_eps,cube],destination=proc)
                            proc_status[proc]=='busy' # Mark processor busy
                            self.cube_grad[cube]='calc %d' %proc # Mark that the cube is being calculated
                            break
            #
            # If we get here there are no free processors
            #
            # Run the next job here and then wait
            #
            for cube in cubes:
                if not self.cube_grad.has_key(cube):
                    score=self.calc_cube(cube,master=True)
                    self.cube_grad[cube]=score
                    to_get=to_get-1
                    break
            #
            # Wait for messages for all jobs
            #
            if self.options.cluster:
                for cube in sorted(self.cube_grad.keys()):
                    if type(self.cube_grad[cube]) is type('text'):
                        node=int(self.cube_grad[cube].split()[-1])
                        score=pypar.receive(node)
                        to_get=to_get-1
                        proc_status[proc]='free'
                        self.cube_grad[cube]=score
            #
            # check if all is done
            #
            txt='Cube scores to get: %4d' %to_get
            print (1+len(txt))*'\b'+txt,
            sys.stdout.flush()
            #
            if len(self.cube_grad.keys())==len(cubes):
                all_done=True
                for cube in sorted(self.cube_grad.keys()):
                    if type(self.cube_grad[cube]) is type('text'):
                        all_done=False
        return self.cube_grad

    #
    # -----
    #
    
    def get_error(self,titgroup,residue,atom,now_array):
        """Get the error on the prediction of a specific ghost
        This function handles all ranges, <, > etc."""
        exp_value=self.exp_ghosts[titgroup][residue][atom]
        calc_value=now_array[titgroup][residue][atom]
        experror=self.exp_ghosts[titgroup][residue][atom+'_error']
        error,satisfied,abs_sat,tot_restraints,tot_abs,real_error=self.get_error_sub(exp_value,calc_value,experror,atom)
        return error,satisfied,abs_sat,tot_restraints,tot_abs
        
    #
    # ---------- 
    #
        
    def get_error_sub(self,exp_value,calc_value,experror,atom,calc_error=0):
        """
        # Subroutine for calculating the error of a given calculated value
        We have a calculated value with an error and an experimental value with an error
        Occasionally we have a range for the experimental value
        """
        error=0.0
        real_error=0.0
        #
        # Internal counters
        #
        tot_restraints=0
        satisfied=0
        tot_abs=0
        abs_sat=0

        #
        if exp_value[0]=='q':
            exp_value=exp_value[1:]
        exp_shift=0.0
        #
        if len(exp_value.split(';'))==2:
            #
            # Deal with errors if the experimental value is a range
            #
            within=True
            totalval=0.0
            errors=[]
            for term in exp_value.split(';'):
                if term[0]=='<':
                    if (calc_value-calc_error)<float(term[1:]):
                        pass
                    else:
                        within=False
                        errors.append((calc_value-calc_error)-float(term[1:]))
                elif term[0]=='>':
                    if (calc_value+calc_error)>float(term[1:]):
                        pass
                    else:
                        within=False
                        errors.append((calc_value+calc_error)-float(term[1:]))
                totalval=totalval+float(term[1:])
            #
            # Were any of the limits violated?
            #
            if within:
                satisfied=satisfied+1
                exp_shift=float(totalval/2.0)
            else:
                real_error=sum(errors)/float(len(errors))
                error=abs(min(errors)) # Not correct
            tot_restraints=tot_restraints+1
            #
        elif exp_value[0]=='<':
            #
            # Less than
            #
            if (calc_value-calc_error)<exp_value:
                satisfied=satisfied+1
                exp_shift=float(exp_value[1:])
            else:
                real_error=(calc_value-calc_error)-exp_value
                error=abs(real_erorr)
            tot_restraints=tot_restraints+1
            
        elif exp_value[0]=='>':
            #
            # greater than
            #
            if calc_value>exp_value:
                satisfied=satisfied+1
                exp_shift=float(exp_value[1:])
            else:
                error=abs(calc_value-exp_value)
            tot_restraints=tot_restraints+1
        elif exp_value=='absent':
            if self.options.use_absent:
                #
                # No ghost
                #
                if (atom=='H' and (calc_value-calc_error)<0.03) or (atom=='N' and (calc_value-calc_error)<0.1):
                    satisfied=satisfied+1
                    abs_sat=abs_sat+1
                else:
                    if atom=='N':
                        real_error=(calc_value-calc_error)-0.1
                        error=abs(real_error)
                    else:
                        real_error=(calc_value-calc_error)-0.03
                        error=abs(real_error)
                #
                tot_restraints=tot_restraints+1
                tot_abs=tot_abs+1
        else:
            #
            # Normal value
            #
            diff=abs(float(exp_value)-calc_value)
            if diff<(experror+calc_error):
                satisfied=satisfied+1
                exp_shift=float(exp_value)
            else:
                error=abs(diff-experror-calc_error)
                real_error=float(exp_value)-calc_value-experror-calc_error
            tot_restraints=tot_restraints+1
        return error,satisfied,abs_sat,tot_restraints,tot_abs,real_error
    
    #
    # ------
    #
    
    def get_score(self,epsmap=None):
        """Do the PBE calculation and compare with the experimental values"""
        if epsmap is None:
            epsmap=self.epsmap
        #
        # Set the epsmap
        #
        if self.options.epsmap:
            for cube in self.cube_eps.keys():
                epsmap.set_cubeeps(cube,self.cube_eps[cube])
        #
        # Get the ghosts
        #
        calc_ghosts=self.get_spans(epsmap=epsmap)
        self.current_ghosts=calc_ghosts.copy()
        #
        # The score is sum of all unsatisfied restraints - score should be minimized
        score=0.0
        sum_exp=0.0
        satisfied=0
        tot_restraints=0
        abs_sat=0
        tot_abs=0
        #
        # Calculate the error
        #
        self.errors={}
        for titgroup in sorted(self.exp_ghosts.keys()):
            self.errors[titgroup]={}
            for residue in sorted(self.exp_ghosts[titgroup].keys()):
                for atom in ['H','N']:
                    if self.exp_ghosts[titgroup][residue].has_key(atom):
                        #
                        error,x_satisfied,x_abs_sat,x_tot_restraints,x_tot_abs=self.get_error(titgroup,residue,atom,calc_ghosts)
                        satisfied=satisfied+x_satisfied
                        abs_sat=abs_sat+x_abs_sat
                        tot_restraints=tot_restraints+x_tot_restraints
                        tot_abs=tot_abs+x_tot_abs
                        #
                        pred=calc_ghosts[titgroup][residue][atom]
                        if self.options.verbose:
                            expshift=self.exp_ghosts[titgroup][residue][atom]
                            print '%13s %13s %2s pred: %6.3f exp: %s %5.3f' %(titgroup, residue, atom, pred, expshift,score)
        #print 'Sum of experimental shifts: %6.3f' %sum_exp
        txt='Restraints satisfied: %5.2f%% (%4d), total # of restraints: %4d. Non-absent restraints satisfied %5.2f%% (%4d/%4d)' %(float(satisfied)/tot_restraints*100.0,satisfied,tot_restraints, float(satisfied-abs_sat)/(tot_restraints-tot_abs)*100.0,satisfied-abs_sat,tot_restraints-tot_abs)
        if self.options.verbose:
            if not epsmap:
                print 'Dielectric constant: %5.2f, %s' %(self.eps,txt)
            else:
                print txt
        else:
            print self.eps,
            import sys
            sys.stdout.flush()
        return 100.0-float(satisfied)/tot_restraints*100.0
    
    #
    # ----
    #
    def get_charge_center(self,group):
        """Get the center of the charge for a group"""
        import string
        if group.split('-')>1:
            titgroup=string.strip(group.split('-')[0])
        else:
            titgroup=group
        tit_type=titgroup.split(':')[-1]
        residue=self.PI.resnum(titgroup)
        centers=self.PI.calculate_titgroup_center(residue)
        if centers.has_key(residue):
            if centers[residue].has_key(tit_type):
                return centers[residue][tit_type]
        raise Exception,'Could not calculate titgroup center'
        
        
    #
    # -----
    #
    
    def get_starting_epsmaps(self):
        """Get the starting dielectric maps"""
        print 'Calculating starting epsmap'
        Emap_ch=False
        if self.options.epsmap:
            self.options.epsmap=False
            Emap_ch=True
        thisdir=self.run_APBS(titgroup=self.titgroups[0],epsmap=None,APBSauto=False,deletedir=False)
        import os, shutil
        for fn in ['xdiel.dx','ydiel.dx','zdiel.dx']:
            shutil.move(os.path.join(thisdir,fn),os.path.join(self.topdir,'start_values',fn))
        # Delete the APBS tempdir
        import shutil
        shutil.rmtree(thisdir)
        print 'Initial dielectric maps generated'
        if Emap_ch:
            # Make sure that epsmap is set to true again.
            self.options.epsmap=True
        return
        
    #
    # ------
    #
        
    def get_spans(self,epsmap=None):
        """Get the spans for the titgroups we are benchmarking"""
        ghosts={}
        for titgroup in sorted(self.titgroups):
            ghosts[titgroup]={}
            #
            grid=None
            if self.options.pbe and not self.options.focus:
                grid=self.run_APBS(titgroup=titgroup,epsmap=epsmap)
            #
            if self.options.verbose:
                print 'TITRATABLE GROUP',titgroup
                print 'Residue  CS Nitrogen    CS Hydrogen    CS H-alpha'
            #
            # Below is for focussing
            #
            for residue in sorted(self.exp_ghosts[titgroup].keys()):
                # Run a calculation for each residue - focussing very closely on the bond
                if self.options.pbe and self.options.focus:
                    grid=self.run_APBS(titgroup=titgroup,epsmap=epsmap,focus=self.PI.GetPosition(residue+':N'))
                #
                if self.exp_ghosts[titgroup][residue].has_key('N'):
                    #
                    # N
                    #
                    atomname=residue+':N'
                    dCS_N=self.get_Ghost(atom=atomname,grid=grid,residue=residue,titgroup=titgroup)
                    ghosts[titgroup][residue]={'N':dCS_N}
                #
                if self.exp_ghosts[titgroup][residue].has_key('H'):
                    #
                    # H
                    #
                    atomname=residue+':H'
                    dCS_H=self.get_Ghost(atom=atomname,grid=grid,residue=residue,titgroup=titgroup)
                    if not ghosts[titgroup].has_key(residue):
                        ghosts[titgroup][residue]={}
                    ghosts[titgroup][residue]['H']=dCS_H
                    #
                if self.exp_ghosts[titgroup][residue].has_key('HA'):
                    #
                    # HA
                    #
                    atomname=residue+':HA'
                    dCS_HA=self.get_Ghost(atom=atomname,grid=grid,residue=residue,titgroup=titgroup)
                    if not ghosts[titgroup].has_key(residue):
                        ghosts[titgroup][residue]={}
                    ghosts[titgroup][residue]['HA']=dCS_H
                #
                # Delete vgrid when focussing
                #
                if self.options.pbe and self.options.focus:
                    delete_vgrid(grid)
            #
            # Delete vgrid when not focussing
            #
            if self.options.pbe and not self.options.focus:
                delete_vgrid(grid)
        return ghosts
        
    #
    # -----
    #    
    
    def get_interaction(self,other_titgroup,grid=None,epsmap=None,residue=None,titgroup=None):
        """Calculate the charge-charge interaction between other_titgroup and titgroup"""
        
        return
    
    
    def get_Ghost(self,atom=None,grid=None,epsmap=None,residue=None,titgroup=None):
        """Calculate Ghosts"""
        if self.options.pbe:
            potdiff=self.get_potdiff(atom,grid)
            #factors=0.025692*ppm_au*1E6/1E-10*NSP[atom.split(':')[-1]]
            factors=0.00049963232399999997*NSP[atom.split(':')[-1]]
            dCS=potdiff*factors
        else:
            import pKa.pKD_tools
            titpos=self.get_charge_center(titgroup)
            titgroup_type=pKa.pKD_tools.get_titgroup_type_from_titgroup(titgroup)
            charge=self.PI.titgroups[titgroup_type]['charge']
            dCS=self.get_span_Coulomb(residue,titpos,charge,self.eps,atom)
        return dCS
        
    #
    # -----
    #
        
    def get_potdiff(self,atom,grid):
        """Get the potential difference between two atoms for ghost titrations"""
        if atom.split(':')[-1]=='N':
            atom2=self.PI.PreviousResidue(self.PI.resid(atom))+':C'
        elif atom.split(':')[-1]=='H':
            atom2=atom[:-1]+'N'
        elif atom.split(':')[-1]=='HA':
            atom2=atom[:-2]+'CA'
        else:
            raise Exception('Unknown atom: %s' %atom)
        #
        potdiff=-(self.get_potential(atom,grid)-self.get_potential(atom2,grid))
        return potdiff/(self.PI.dist(atom,atom2))
        
    #
    # -----
    #
        
    def get_potential(self,atom,grid):
        """Get the potential at this atom"""
        position=self.PI.GetPosition(atom)
        return self.get_potential_atpos(position,grid)
    #
    # -----
    #
              
    def get_potential_atpos(self,position,grid):
        """Get the potential at a specific position"""
        inval=0.0
        pos=[position[0],position[1],position[2]]
        ret,potential=Vgrid_value(grid,pos,inval)
        return potential
        
        
    #
    # -------
    #
        
    def run_APBS(self,titgroup=None,charge=None,epsmap=None,focus=None,APBSauto=True,deletedir=True):
        """Run APBS in mg-auto mode to get a map focussed on the protein"""
        #
        # Potname
        #
        import tempfile
        thisdir=tempfile.mkdtemp(dir='/tmp')
        import os
        potname=os.path.join(thisdir,'%s_potential' %(titgroup))
        #
        # Get charge
        #
        import pKa.pKD_tools
        titgroup_type=pKa.pKD_tools.get_titgroup_type_from_titgroup(titgroup)
        charge=self.PI.titgroups[titgroup_type]['charge']
        #
        # Get the PQR file
        #
        import Protool.assign_parameters
        self.PQR=Protool.assign_parameters.assign_parameters(self.PI)
        self.PQR.clear_all_charges_and_radii()
        self.PQR.set_all_radii()
        titpos=self.get_charge_center(titgroup)
        self.PI.add_atom(uniqueid='Z:999:CHA',atomnumber=99999,atomname='CHA',    
                            chainid='Z',residuename='DUM',residuenumber='999',
                            xcoord=titpos[0],ycoord=titpos[1],zcoord=titpos[2],update=1,BFACTOR=None,OCCUPANCY=None,CHARGE=charge,RADIUS=1.55,tag=None)
        pqrfile=os.path.join(thisdir,'currentPQR.pqr')
        self.PI.writepqr(pqrfile)
        self.PI.remove_atom('Z:999:CHA')
        #
        # Get the coarse and fine centers
        #
        self.coarsecent,extent=self.PI.get_center_of_coords()
        self.coarsedim=1.5*extent
        if focus is None:
            self.finedim=extent
            self.finecent=self.coarsecent
        else:
            self.finecent=[focus[0],focus[1],focus[2]] 
            ngrid=65
            self.finedim=[ngrid*0.2,ngrid*0.2,ngrid*0.2]
        #
        # Run APBS
        #
        import os
        os.chdir(thisdir)
        infile=os.path.join(thisdir,'apbs.in')
        fp = open(infile, "w")
        text = "read\n"
        text += "    mol pqr %s\n" %pqrfile
        #
        #
        #
        if self.options.epsmap:
            name='usemap'
            dirname=thisdir
            names=self.epsmap.set_and_write_epsmap(dirname,name)
            text += "    diel dx %s %s %s\n" %(names[0],names[1],names[2])
            #text += "    kappa dx kappa%d.dx\n"%i
            #text += '    diel dx '
            #for n in ['x','y','z']:
            #    text=text+" %s " %(os.path.join(self.topdir,'start_values','%sdiel.dx' %n))
            #text=text+" \n"
            
        text += "end\n"
        #
        # Start of ELEC section
        #
        text += "elec\n"
        if APBSauto:
            text += "    mg-auto\n"
        else:
            text=text+"    mg-manual"
            
        text += "    dime 65 65 65\n"
        #text += "    grid %.2f %.2f %.2f\n"%(dim, dim, dim)
        #text += "    gcent %.3f %.3f %.3f\n"%(self.center[0], self.center[1], self.center[2])
       

        if APBSauto:
            text += "    cglen %.4f %.4f %.4f\n" % (self.coarsedim[0], self.coarsedim[1], self.coarsedim[2])
            text += "    cgcent %.3f %.3f %.3f\n" %(self.coarsecent[0],self.coarsecent[1],self.coarsecent[2])
            text += "    fglen %.4f %.4f %.4f\n" % (self.finedim[0], self.finedim[1], self.finedim[2])
            text += "    fgcent %.3f %.3f %.3f\n" %(self.finecent[0],self.finecent[1],self.finecent[2])
        else:
            text += "    glen %.4f %.4f %.4f\n" % (self.coarsedim[0], self.coarsedim[1], self.coarsedim[2])
            text += "    gcent %.3f %.3f %.3f\n" %(self.coarsecent[0],self.coarsecent[1],self.coarsecent[2])
        #
        text += "    mol 1\n"
        text += "    lpbe\n"
        text += "    ion charge 1 conc 0.00 radius 2.0\n"
        text += "    ion charge -1 conc 0.00 radius 2.0\n"
        if self.options.smoothe:
            text += "    srfm smol\n" # smoothe the surface
        else:
            text += "    srfm mol\n"  # do not smoothe the dielectric boundary!
        text += "    chgm spl0\n"
        text += "    srad 1.4\n"
        text += "    swin 0.3\n"
        text += "    sdens 10.0\n"
        text += "    temp 298.15\n"
        #text += "    usemap kappa 1\n"
        text += "    pdie %5.2f\n" %self.eps
        text += "    sdie %5.2f\n" %self.options.wateps
        if self.options.epsmap:
            text += "    usemap diel 1-3\n"
        text += "    bcfl sdh"
        text += "    calcenergy no\n"
        text += "    calcforce no\n"
        text += "    write dielx dx xdiel\n"
        text += "    write diely dx ydiel\n"
        text += "    write dielz dx zdiel\n"
        text += "    write pot dx %s\n" %potname
        #text += "    write kappa dx kappa\n"

        text += "end\n"
            
        fp.write(text)
        fp.close()
        # run APBS #
        if not options.verbose:
            import subprocess
            outputStream = open('/dev/null', 'w+')
            process = subprocess.Popen(['%s %s' %(APBSpath,infile)], stdout=outputStream, stderr=subprocess.STDOUT, shell=True)
            process.wait()
            outputStream.close()
            os.system("rm io.mc")
        else:
            import os
            os.system('%s %s' %(APBSpath,infile))
        #
        data=[]
        value=0.0
        startVio()
        grid = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, data)
        Vgrid_readDX(grid, "FILE", "ASC", "", potname+'.dx')
        if deletedir:
            #
            # Delete the dir
            #
            import shutil
            shutil.rmtree(thisdir)
            return grid
        else:
            return thisdir
    
    #
    # -----
    #
        
    def get_angle(self,residue,titpos,atom_type):
        """Get the cos(angle)and distance between a bond and a charge"""
        if atom_type=='N':
            try:
                prev_res=self.PI.PreviousResidue(residue)
            except:
                return None,None
            atoms=[residue+':N',prev_res+':C']
        elif atom_type=='H':
            atoms=[residue+':H',residue+':N']
        elif atom_type=='HA':
            atoms=[residue+':HA',residue+':CA']
        else:
            print 'Atom type %s unknown' %str(atom_type)
            raise Exception
        #
        # Check that both atoms are present
        #
        if not self.PI.atoms.has_key(atoms[0]) or not self.PI.atoms.has_key(atoms[1]):
            raise Exception('Atoms are not present: %s %s' %(atoms[0],atoms[1]))
        #
        # Get the bond vector
        #
        bond_vector=self.PI.GetPosition(atoms[1])-self.PI.GetPosition(atoms[0])
        #
        # Vector to charge
        #
        charge_vector=titpos-self.PI.GetPosition(atoms[0])
        #
        # Get angle
        #
        import numpy
        dp=numpy.dot(bond_vector,charge_vector)
        cos_angle=dp/(length(bond_vector)*length(charge_vector))
        #
        # Get the distance
        #
        dist=length(charge_vector)
        return cos_angle,dist

    #
    # ----
    #
        
    def get_span_Coulomb(self,residue,titpos,charge,deff,atom='N'):
        """When given a residue calculate the backbone and possible side chain polarisation due to
        the titratable charge at position titpos. The titratable group has a charge = charge and
        interacts with the residue with an effective dielectric constant deff        """
        import numpy
        import os
        atom_type=atom.split(':')[-1]
        cos_angle,dist=self.get_angle(residue=residue,titpos=titpos,atom_type=atom_type)
        if not cos_angle or not dist:
            return None,None
        #
        # Calculate electric field
        #
        import math
        e=1.602E-19
        e0=8.85E-12
        CS=1.0/(4*math.pi*e0*deff)*e/(dist*1E-10)**2*ppm_au*NSP[atom_type]*1E6*cos_angle*charge
        return CS
        
    
#
# ------------
#

def build_hydrogens(pdbfile):
    """Build hydrogens with Protool"""
    import Protool
    PI=Protool.structureIO()
    PI.readpdb(pdbfile)
    PI.RemoveALT()
    #
    #
    # Check if His 15 is called HID
    #
    #print sorted(PI.residues.keys())
    for residue in sorted(PI.residues.keys()):
        if PI.resname(residue)=='HID':
            for atom in PI.residues[residue]:
                PI.atoms[atom]['RESNAME']='HIS'
        elif PI.resname(residue)=='CYX':
            for atom in PI.residues[residue]:
                PI.atoms[atom]['RESNAME']='CYS'
    PI.Update()
    PI.Remove_All_NonAminoAcids()
    PI.remove_all_hydrogens()
    #
    # Keep only the first chain
    #
    chains=sorted(PI.chains.keys())
    if len(chains)>1:
        for chain in chains[1:]:
            print 'Deleting chain %s' %chain
            PI.remove_chain(chain,update=True)
    #
    PI.setcha(PI.residues.keys(),'') # Change chain ID
    #
    import Protool.hydrogens
    H=Protool.hydrogens.Hydrogens(PI)
    H.build_all_HNs()
    H.build_all_HAs()
    outfile=pdbfile[:-4]+'_H.pdb'
    PI.writepdb(outfile)
    return outfile

#
# -----
#

def dielscan(options,args):
    """Do a scan of the dielectric constant"""
    indi_errors={}
    results={}
    all_ghosts={}
    for method in ['coulomb','PBE']:
        all_ghosts[method]={}
        txt='Method: %s' %method
        print method,
        import sys
        sys.stdout.flush()
        for absent in ['use']: #,'dontuse']:
            options.smoothe=False
            if method=='PBE':
                options.pbe=True
            elif method=='sPBE':
                options.pbe=True
                options.smoothe=True
            else:
                options.pbe=False
            if absent=='use':
                options.use_absent=True
            else:
                options.use_absent=False
            xs=[]
            ys=[]
            best_eps=None
            best_score=10000.0
            for eps in range(5,300,5): 
                xs.append(eps/10.0)
                X=iterative_epsmap(eps/10.0,options)
                score=X.get_score()
                if score<best_score:
                    best_eps=eps
                    best_score=score
                ys.append(score)
                #
                # Store all the calculated ghosts
                #
                import copy
                all_ghosts[method][eps]=copy.deepcopy(X.current_ghosts)
    #
    # Save ghosts
    #
    import os
    use_name=os.path.split(options.pdb)[1]
    picklefile=os.path.join(topdir,'all_ghosts_%s.pickle' %use_name)
    fd=open(picklefile,'w')
    import pickle
    pickle.dump(all_ghosts,fd)
    fd.close()
    return

    #
    # -----
    #
 
def main(options,args):
    """Main program"""
    #
    # Make the start_values dir
    #
    ID=0
    if options.cluster:
        import pypar
        ID=pypar.rank()
    #
    import os
    startdir=os.path.join(topdir,'start_values')
    if ID==0:
        if not os.path.isdir(startdir):
            os.mkdir(startdir)
    #
    # Adjust the PDB file name
    #
    options.pdbfile=os.path.join(topdir,options.pdb)
    #
    # Make the starting epsmaps
    #
    if ID==0 and options.epsmap:
        X=iterative_epsmap(options.eps,options,init_epsmap=True)
        X.get_starting_epsmaps()
    if options.cluster:
        pypar.barrier()
    #
    # See what we should do
    #
    if options.dielscan:
        # dielectric scan
        count=0
        if options.allPDBs:
            import os
            failed=[]
            pdbs=os.listdir(topdir)
            os.chdir(topdir)
            for pdb in sorted(pdbs):
                if pdb[-4:]=='.pdb' and pdb[-6:]!='_H.pdb':
                    pdbname=pdb
                    full_pdbpath=os.path.join(topdir,pdbname)
                    #
                    # Should we consider this MD snapshot? (if this is a snapshot?)
                    #
                    if options.MD!=0:
                        number=int(pdb.split('.')[-2])
                        if number%options.MD!=0:
                            continue
                    #
                    # Should we build hydrogens?
                    #
                    try:
                        if options.buildHs:
                            print 'building hydrogens for %s' %pdb
                            pdbname=build_hydrogens(full_pdbpath)
                    except:
                        print 'Building hydrogens failed. Skipping PDB'
                        continue
                    options.pdb=pdbname
                    pdb_short=os.path.split(pdbname)[1]
                    #
                    count=count+1
                    txt= 'Doing dielscan for %25s, which is number %4d' %(pdb,count)
                    print txt,
                    import sys, os
                    sys.stdout.flush()
                    resultfile=os.path.join(topdir,'all_ghosts_%s.pickle' %pdb_short)
                    workfile=os.path.join(topdir,'%s_working.status' %pdb_short)
                    if not os.path.isfile(workfile):
                        if not os.path.isfile(resultfile):
                            fd=open(workfile,'w')
                            fd.write('yes, indeed\n')
                            fd.close()
                            #try:
                            dielscan(options,args)
                            print 'Success!'
                            import os
                            os.remove(workfile)
                            #except:
                            #    failed.append(pdb)
                            #    print 'Failed :-('
                            #    #raise Exception()
                        else:
                            print 'Already done!'
                    else:
                        print 'Being processed'
            #
            print 'Failed for %3d pdb files' %(len(failed))
            print sorted(failed)

        else:
            #
            # Just do a scan for a single PDB
            #
            if options.buildHs:
                print 'building hydrogens for %s' %options.pdb
                pdbname=build_hydrogens(options.pdb)
                options.pdb=pdbname
            dielscan(options,args)
        return
    #
    if options.epsmap:
        # Optimising the epsmap
        X=iterative_epsmap(options.eps,options)
        if ID==0:
            if options.cubescan:
                X.cube_scan()
            elif options.LM:
                X.LM_find_epsmap()
            elif options.fillact>0.0:
                X.fill_actsite(options.fillact,options.epsact,80.0)
            else:
                X.find_epsmap()
        else:
            X.slave_wait()
        return
    #
    # Make experimental data
    #
    if options.make_expdata:
        #
        # Calculate ghosts and save them
        #
        import os
        X=iterative_epsmap(options.eps,options)
        X.titgroups=X.PDB_titgroups[:]
        X.exp_ghosts={}
        for tg in X.titgroups:
            X.exp_ghosts[tg]={}
            for residue in sorted(X.PI.residues.keys())[1:]:
                X.exp_ghosts[tg][residue]={'N':0.0}
                if X.PI.resname(residue)!='PRO':
                    X.exp_ghosts[tg][residue]['H']=0.0
        ghosts=X.get_spans()
        #
        # reformat
        #
        wg={}
        for tg in ghosts.keys():
            wg[tg]={}
            for residue in ghosts[tg].keys():
                wg[tg][residue]={}
                for atom in ghosts[tg][residue].keys():
                    import types
                    if not type(ghosts[tg][residue][atom]) is types.TupleType:
                        wg[tg][residue][atom]='%5.3f' %(float(ghosts[tg][residue][atom]))
        #
        os.chdir(topdir)
        print 'Trying to write file to local dir'
        import os
        print os.getcwd()
        print options.restraints
        fd=open(options.restraints,'w')
        import pickle
        pickle.dump(wg,fd)
        fd.close()
        print 'Simulated experimental restraints written to %s' %options.restraints
        return
    #
    # If nothing else specified then we just calculate the score and save the predictions
    #
    X=iterative_epsmap(options.eps,options)
    score=X.get_score()
    import os
    os.chdir(topdir)
    fd=open('current_ghosts.pickle','w')
    import pickle
    pickle.dump(X.current_ghosts,fd)
    fd.close()
    #
    print 'Dielectric constant: %5.2f, score: %5.2f %%' %(options.eps,score)
    return
#
# -----
#
if __name__=='__main__':
    import optparse
    parser=optparse.OptionParser()
    parser.add_option('-e','--eps',dest='eps',type='float',action='store',default=8.0,
        help='Protein dielectric constant (only valid if epsmap is not specified) Default: %default')
    parser.add_option('-w','--wateps',dest='wateps',type='float',action='store',default=80.0,
        help='Solvent dielectric constant Default: %default')
    parser.add_option('-s','--sidelength',dest='sidelength',type='float',action='store',default=5.0,
        help='Side length of cubes in epsmap. Default: %default')
    parser.add_option('-z','--stepsize',dest='stepsize',type='float',action='store',default=1.0,
        help='Stepsize for minimization. Default: %default')
    #
    parser.add_option('--pbe',dest='pbe',action='store_true',default=False,
        help='Use PBE solver. Default: %default')
    parser.add_option('--focus',dest='focus',action='store_true',default=False,
        help='Focus on each peptide bond when calculating the CS. Default: %default')
    parser.add_option('--smoothe',dest='smoothe',action='store_true',default=False,
        help='Smoothe the PBE surface. Default: %default')
    #
    # PDB file
    #
    parser.add_option('--pdb',dest='pdb',default='2LZT_H.pdb',
                      help='PDB file. Default: %default')
    parser.add_option('--allpdbs',dest='allPDBs',default=False,action='store_true',
                      help='Perform action for all PDB files in pwd. Only implemented for dielscan. Default: %default')
    parser.add_option('--buildHs',dest='buildHs',default=False,action='store_true',
                      help='Build HN hydrogens. Default: %default')
    parser.add_option('--MD',dest='MD',default=0,action='store',type='int',help='Integer multiple of which to select MD snapshots. Default: %default')
    parser.add_option('--pdblist',dest='pdblist',default='',
                        help='File containing list of PDB files')
    #
    # Experimental data
    #
    parser.add_option('-x','--expdata',dest='restraints',default='Experimental_ghosts.txt',
                      help='Experimental data restraints. Default: %default')
    parser.add_option('--use_absent',dest='use_absent',action='store_true',default=False,
                      help='Use absent ghosts to calculate score. Default: %default')
    #
    parser.add_option('--dielscan',dest='dielscan',action='store_true',default=False,
        help='Do a scan with different values for the protein dielectric constant. Default: %default')
    parser.add_option('--cubescan',dest='cubescan',action='store_true',default=False,
        help='Do a scan of all cubes seeing which cube gives the biggest improvement when its eps is halved. Default: %default')
    parser.add_option('--fillact',dest='fillact',action='store',type='float',default=0.0,
        help='If different from zero then fill active site with epsact to this depth. Default: %default')
    parser.add_option('--epsact',dest='epsact',action='store',type='float',default=8.0,
        help='dielectric constant to fill active site with. Default: %default')
    #
    # Generate mock experimental data with noise
    #
    parser.add_option('--make_expdata',dest='make_expdata',action='store_true',default=False,
                      help='Construct dictionary with mock experimental data. Default: %default')
    #
    # Optimization options
    #
    parser.add_option('--epsmap',dest='epsmap',action='store_true',default=False,
        help='Use 3D epsmap in PBE calculations. Default: %default')
    parser.add_option('--opteps',dest='opteps',action='store_true',default=False,
        help='Optimze the 3D epsmap in PBE calculations. Default: %default')
        
    parser.add_option('--noLM',dest='LM',action='store_false',default=True,
        help='Use LM minimization. Default: %default')
    #
    parser.add_option('--starteps',dest='starteps',type='string',default='',
                      help='Starting epsmap (in pickle format)')
    #
    # Execution options
    #
    parser.add_option('-v','--verbose',dest='verbose',action='store_true',default=False,
        help='Verbose. Provide more output')
    parser.add_option('--par',dest='cluster',action='store_true',default=False,
                      help='Execute code in parallel with pypar. Default: %default')
    #
    # Parse the input
    #
    (options, args,) = parser.parse_args()
    #
    #
    if (options.cubescan or options.opteps) and not options.epsmap:
        parser.error('You must specify --epsmap with --cubescan and --opteps')
    #
    if options.epsmap and not options.pbe:
        parser.error('You must specify --pbe with --epsmap')
    
    main(options,args)
    import sys
    sys.exit(0)
    
    

                
