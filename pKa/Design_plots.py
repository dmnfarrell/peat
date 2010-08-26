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

import math
basedir='/home/nielsen/pKa-design/done_distnummuts'
two_dir='/home/nielsen/pKa-design/design_two'

def get_prefix(name):
    #
    # Set the prefix
    #
    if name.find('distance_nummuts__')!=-1:
        prefix=name.split('distance_nummuts__')[-1]
        prefix=prefix.split('_')[0]
    else:
        prefix=name.split('designtwo__')[-1]
        prefix=prefix.split('_')[0]
    return prefix

#
# -----
#

def average(l):
    sum=0.0
    for i in l:
        sum=sum+i
    #print 'Sum: %f, dividing by %d' %(sum,len(l))
    avg=sum/float(len(l))
    #
    # Variance
    #
    N=len(l)
    sum2=0.0
    for i in l:
        sum2=sum2+math.pow((i-avg),2)
    if N>1:
        variance=sum2/float(N-1)
    else:
        variance=0
        print 'Warning - N was 1'
    std_dev=math.sqrt(variance)
    return avg,variance,std_dev

#
# ----
#

def read_actsit_def(filename='/home/nielsen/pKa-design/done_distnummuts/active_sites'):
    """Read the active site definition"""
    d={}
    fd=open(filename)
    line=fd.readline()
    current_name=None
    while line:
        if line[0]=='*':
            current_name=line[1:].strip()
            d[current_name]=[]
        else:
            if current_name:
                d[current_name].append(line.strip())
        line=fd.readline()
    fd.close()
    return d


#
# ----
#

def get_sum(l):
    s=0.0
    for i in l:
        s=s+i
    return s

#
# ----
#

def get_files(only=None,type=None,selection=None):
    #
    # Get the names of all the files
    #
    loadfiles=[]
    import os
    if type=='two':
        dir=two_dir
    else:
        dir=basedir
    dirs=os.listdir(dir)
    #
    # Loop over all subdirs
    #
    for subdir in dirs:
        if only:
            if not subdir in only:
                continue
        if os.path.isdir(os.path.join(dir,subdir)):
            rdir=os.path.join(dir,subdir)
            #
            # Find the dist_nummuts _MC files
            #
            files=os.listdir(rdir)
            import string
            for file in files:
                if type=='two':
                    if file.find(selection)!=-1:
                        loadfiles.append(os.path.join(rdir,file))
                else:
                    if string.find(file,'pdb_%s' %type)!=-1:
                        loadfiles.append(os.path.join(rdir,file))
    return loadfiles

#
# ------------------
#

def main():
    #
    # Do all the analyses we want
    #
    print
    print 'Design_plots.py: Do all analyses of the Design_pKa runs'
    print
    print 'Usage Design_plots.py [files] <type>'
    print
    #
    # dpKa vs. distance from active site & number of mutations
    #
    import sys
    #
    # Get the type
    #
    type=sys.argv[-1]
    if type=='two':
        #
        # Analysing a single group
        #
        files=get_files(sys.argv[1:-2],type,sys.argv[-2])
    else:
        #
        # Get the files
        #
        files=get_files(sys.argv[1:-1],type)
    #
    # If not files then exit
    #
    if files==[]:
        print 'Error: Did not find any files to match criteria'
        return
    #
    # Prepare the data matrix
    #
    raw_data={}
    max_dist=25
    max_muts=20
    distance_range=range(max_dist+1)
    nummuts_range=range(1,max_muts+1)
    for num in nummuts_range:
        raw_data[num]={}
        for dist in distance_range:
            raw_data[num][dist]=[0.0]
    
    #
    # Loop over all the files
    #
    added=0
    big_dict={}
    tot_target={}
    for file in files:
        if file[-5:]=='.lock':
            continue
        print 'Processing %s' %file
        try:
            import pickle
            fd=open(file)
            d=pickle.load(fd)
            fd.close()
        except:
            continue
        #
        # Set the prefix
        #
        prefix=get_prefix(file)
        #
        # -----------------------------------
        #
        # Loop over all the design-data
        #
        targets=d.keys()
        targets.sort()
        for target in targets:
            #if target!=':0231:ASP':
            #    continue
            #
            # pdbfile and wt_full are not interesting
            #
            if target=='pdbfile' or target=='wt_full':
                continue
            target_pka=d[target]
            designs=target_pka.keys()
            designs.sort()
            if designs==['pKa out of range']:
                continue
            #
            # Loop over each design (normally +20 and -20 for Design_dist_nummuts)
            #
            for design in designs:
                #if design!='m20':
                #    continue
                try:
                    nummuts=target_pka[design].keys()
                except:
                    #print 'Skipping:',target_pka[design]
                    continue
                nummuts.sort()
                for num in nummuts:
                    dist_cutoffs=target_pka[design][num].keys()
                    for cutoff in dist_cutoffs:
                        #text='%15s %4s #muts: %2d, dist_co: %5.2f, sols:' %(target,design,num,float(cutoff))
                        #print text
                        #
                        # Make sure we have a bin for the this distance cutoff
                        #
                        #if not raw_data[num].has_key(cutoff):
                        #    raw_data[num][cutoff]=[]
                        #
                        # Loop over all solutions and store the dpKa values
                        #
                        sol_dict=target_pka[design][num][cutoff]
                        solutions=sol_dict.keys()
                        #
                        # Loop over all the solutions
                        #
                        for sol in solutions:
                            if sol_dict[sol].has_key(type):
                                dpka=sol_dict[sol][type][target]
                                mutations=sol_dict[sol]['mutations']
                                #
                                # Count the number of mutations
                                #
                                nums=0
                                for mut in mutations:
                                    if mut:
                                        nums=nums+1
                                #
                                # Add the data to the array
                                #
                                # We skip all data points outside the range specified
                                # by max_muts and max_dist
                                #
                                skip=None
                                if not raw_data.has_key(nums):
                                    skip=1
                                if not skip:
                                    if not raw_data[nums].has_key(cutoff):
                                        skip=1
                                if not skip:
                                    raw_data[nums][cutoff].append(dpka)
                                #
                                # Add to the big dictionary
                                #
                                import os
                                
                                tname=prefix+target
                                if not big_dict.has_key(tname):
                                    big_dict[tname]=[]
                                clean_muts=[]
                                for mut in mutations:
                                    if mut:
                                        clean_muts.append(mut)
                                big_dict[tname].append([clean_muts,dpka])

                                #
                                # Keep track of how many we add
                                #
                                added=added+1
                                #print 'Adding: nummuts: %2d, cutoff: %4.1f, dpka: %4.2f' %(nums,cutoff,dpka)
                                #except:
                                #    pass
            #print '--------------------'
    #
    # Read the definition of the active site
    #
    act_site=read_actsit_def()
    #
    # Get properties from the PDB files/wt pKa calculation
    #
    import string, os
    for file in files:
        if file[-5:]=='.lock':
            continue
        prefix=get_prefix(file)
        #
        # Analysis
        #
        print 'Analysing for %s' %prefix
        #
        # Read the PDB file
        #
        
        pdbfile=os.path.join(basedir,prefix[:4],prefix)
        import Protool
        Z=Protool.structureIO()
        Z.readpdb(pdbfile)
        #
        # Get the relative accs
        #
        import WI_tools
        accs=WI_tools.relative_accessibility(pdbfile)
        #
        # Open the wt pKa calc
        #
        import pKaTool.pKaIO as pKaIO
        X=pKaIO.pKaIO(pdbfile)
        pkavals=X.readpka()
        matrix=X.read_matrix()
        for residue in pkavals.keys():
            target=prefix+residue
            if not tot_target.has_key(target):
                tot_target[target]={}
            tot_target[target]['pKa']=pkavals[residue]['pKa']
            elecs=[]
            for other_res in matrix[residue].keys():
                elecs.append(matrix[residue][other_res][0])
            tot_target[target]['elecs']=elecs[:]
            #
            # Insert number of aas
            #
            tot_target[target]['prot_aas']=len(Z.residues.keys())
            #
            # Is this target in the vicinity of the active site?
            #
            tot_target[target]['act_site']=None
            target_res=target.split('pdb')[1]
            target_res=':'+target_res.split(':')[1]
            try:
                target_atoms=Z.residues[target_res]
            except:
                print target_res
                print Z.residues.keys()
                stop
            if act_site.has_key(prefix):
                for act_res in act_site[prefix]:
                    r_act_res=':'+act_res.split(':')[1]
                    for atom2 in Z.residues[r_act_res]:
                        for target_atom in target_atoms:
                            #print 'Comparing',target_atom,atom2
                            if Z.distance(target_atom,atom2)<5.0:
                                tot_target[target]['act_site']='Yes'
            #
            # Insert rel. acc
            #
            if residue[-6:]==':CTERM':
                residue=residue[:-6]
            if residue[-6:]==':NTERM':
                residue=residue[:-6]
            #print accs[residue]['sum']
            tot_target[target]['relacc']=accs[residue]['sum']['rel']
            #print residue,accs[residue]
    print
    print ' All done'
    #
    # How many solutions in total?
    #
    print 'I added %5d solutions to the matrix' %added
    #
    # For each target, what's the maximum dpKa?
    #
    targets=big_dict.keys()
    targets.sort()
    max_dpkas=[]
    all=[]
    actsite_dpkas=[]
    all_actsite_dpkas=[]
    file_dpkas={}
    for target in targets:
        tmp_dpkas=[]
        for solution,dpka in big_dict[target]:
            tmp_dpkas.append(abs(dpka))
            if not file_dpkas.has_key(target[:4]):
                file_dpkas[target[:4]]={}
                file_dpkas[target[:4]]['dpkas']=[]
                file_dpkas[target[:4]]['num_target']=0
                file_dpkas[target[:4]]['max_dpka']=0.0
            #
            # Add the new dpKa
            #
            file_dpkas[target[:4]]['dpkas'].append(abs(dpka))
            
        avg,var,sdev=average(tmp_dpkas)
        print 'Average pKa shift for %25s is %5.2f (%5.2f)' %(target,avg,sdev)
        tmp_dpkas.sort()
        max_dpka=tmp_dpkas[-1]
        max_dpkas.append(max_dpka)
        all=all+tmp_dpkas
        #
        # Store the average and max dpka for each target
        #
        tot_target[target]['avg_dpka']=avg
        tot_target[target]['max_dpka']=max_dpka
        #
        # Set the aa size
        #
        file_dpkas[target[:4]]['prot_aas']=tot_target[target]['prot_aas']
        #
        # Increment the number of targets designed for this protein
        #
        file_dpkas[target[:4]]['num_target']=file_dpkas[target[:4]]['num_target']+1
        #
        # Is is an active site target?
        #
        if tot_target[target]['act_site']:
            actsite_dpkas.append(max_dpka)
            all_actsite_dpkas=all_actsite_dpkas+tmp_dpkas
    #
    # Write the PDB files
    #
    for file in files:
        tf_max=[]
        if file[-5:]=='.lock':
            continue
        prefix=get_prefix(file)
        #
        # Writing Yasara script
        #
        print 'Writing Yasara script for %s' %prefix
        #
        # Read the PDB file
        #
        fd=open('yasara.mcr','w')
        pdbfile=os.path.join(basedir,prefix[:4],prefix)
        fd.write('LoadPDB %s\n' %pdbfile)
        fd.write('ColorAll 606060\n')
        fd.write('Style Stick\n')
        fd.write('HUD Off\n')
        fd.write('HideRes Hoh\n')
        import Protool
        Z=Protool.structureIO()
        Z.readpdb(pdbfile)
        #
        # Zero all B-factors
        #
        #for residue in Z.residues.keys():
        #    for atom in Z.residues[residue]:
        #        Z.atoms[atom]['B-factor']=0.0
        #
        # Loop over all targets and set the colour
        #
        colors={1.0:'Blue',
                2.0:'Cyan',
                3.0:'Green',
                4.0:'Yellow',
                5.0:'Red'}
        for target in tot_target.keys():
            #
            # Collect stats on the max abs(dpka)
            #
            pos=target.find(prefix)
            if pos!=-1:
                if tot_target[target].has_key('max_dpka'):
                    tf_max.append(abs(tot_target[target]['max_dpka']))
            #
            # Write the PDB file
            #
            if pos!=-1:
                resnum=target[pos+len(prefix):]
                resnum=':'+resnum.split(':')[1]
                if Z.residues.has_key(resnum):
                    if tot_target[target].has_key('max_dpka'):
                        col_cutoff=colors.keys()
                        col_cutoff.sort()
                        co=0.5
                        for col in col_cutoff:
                            co=col
                            if tot_target[target]['max_dpka']<col:
                                break
                        colour=colors[co]
                        fd.write('ColorRes %d,%s\n' %(int(resnum[1:]),colour))
                    else:
                        fd.write('ColorRes %d,%s\n' %(int(resnum[1:]),'aaaaaa'))
                else:
                    raise 'Residue not found',target
        #Z.writepdb('BF_dpKa')
        print 'Number of max_pkas in %s is %d' %(prefix,len(tf_max))
        avg,var,sdev=average(tf_max)
        print '%s, average max dpKa %5.1f, sdev: %5.1f' %(prefix,avg,sdev)
    #fd.write('exit\n')
    fd.close()
    #
    # Print all the stats
    #
    print
    print 'Number of targets designed is       : %4d  ' %(len(targets))
    all_targets=len(tot_target.keys())
    print 'Number of targets in total:         : %4d  ' %all_targets
    print '%% designed                         : %5.2f' %(float(len(targets))/float(all_targets)*100.0)
    print
    print 'Number of active site targets       : %5.2f' %(len(actsite_dpkas))
    #
    # Get average Delta pKas
    #
    avg,var,sdev=average(all)
    print 'Average dpKa for all targets is      : %5.2f (%5.2f)' %(avg,sdev)
    avg,var,sdev=average(max_dpkas)
    print 'Average MAX dpKa for all targets is  :%5.2f  (%5.2f)' %(avg,sdev)
    # Max dpka for active sites
    avg,var,sdev=average(actsite_dpkas)
    print 'Average MAX dpKa for active site targets: %5.2f (%5.2f)' %(avg,sdev)
    #
    avg,var,sdev=average(all_actsite_dpkas)
    print 'Average dpKa for actsit target   :%5.2f  (%5.2f)' %(avg,sdev)
    print
    print 'Average dpKa per protein'
    prots=file_dpkas.keys()
    prots.sort()
    for prot in prots:
        avg,var,sdev=average(file_dpkas[prot]['dpkas'])
        num_target=file_dpkas[prot]['num_target']
        aa_size=file_dpkas[prot]['prot_aas']
        num_sol=len(file_dpkas[prot]['dpkas'])
        print 'Average dpKa for %s is               : %5.2f (%5.2f) [#targets %4d, #aas %4d, #sols/target %5.2f]' %(prot,avg,sdev,num_target,aa_size,float(num_sol)/float(num_target))
    #
    # Stats on the types of targets designed
    #
    designed={}
    import pKarun
    Y=pKarun.pKanalyse()
    for target in big_dict.keys():
        rtype=Y.get_residue_type(target)
        if not designed.has_key(rtype):
            designed[rtype]=0
        designed[rtype]=designed[rtype]+1
    des=designed.keys()
    #
    # Look at the targets not designed
    #
    not_designed={}
    all_targets=tot_target.keys()
    all_targets.sort()
    import pKarun
    Y=pKarun.pKanalyse()
    for target in all_targets:
        if not big_dict.has_key(target):
            rtype=Y.get_residue_type(target)
            if not not_designed.has_key(rtype):
                not_designed[rtype]=0
            not_designed[rtype]=not_designed[rtype]+1
    #
    # Stats
    #
    print
    print 'Stats on types of groups designed'
    types=['ASP','GLU','TYR','CYS','CTERM','NTERM','LYS','ARG','HIS']
    types.sort()
    for rtyp in types:
        if designed.has_key(rtyp):
            des=designed[rtyp]
        else:
            des=0
        if not_designed.has_key(rtyp):
            ndes=not_designed[rtyp]
        else:
            ndes=0
        tot=ndes+des
        if tot>0:
            avg='%5.2f' %(float(des)/float(tot)*100.0)
        else:
            avg='NA'
        print '%8s des: %3d notD: %3d, tot: %3d %% designed: %s' %(rtyp,des,ndes,tot,avg)
    #
    # Relation between average dpKa obtained and accessibility, type and electrostatic interactions.
    #
    print
    #
    # Plot of avg dpKa vs. sum of abs electrostatic interactions
    #
    avg_dpka=[]
    max_dpka=[]
    sum_elec=[]
    acc=[]
    for target in all_targets:
        dpkas=[]
        if big_dict.has_key(target):
            for mutants,dpka in big_dict[target]:
                dpkas.append(abs(dpka))
            e_sum=[]
            for elec in tot_target[target]['elecs']:
                e_sum.append(elec)
            #
            max_dpka.append(max(dpkas))
            #
            avg,var,sdev=average(dpkas)
            avg_dpka.append(avg)
            #
            avg,var,sdev=average(e_sum)
            sum_elec.append(get_sum(e_sum))
            #
            # Accessibility
            #
            acc.append(tot_target[target]['relacc'])
        else:
            #print 'No design for',target
            pass
    import dislin_driver
    file=dislin_driver.graf_mult2(acc,[avg_dpka,max_dpka],
                                  title='Effect of solvent exposure',
                                  x_legend='Relative accessibility of target',
                                  y_legend='abs(dpKa)',
                                  legends=['Avg. dpKa','Max. dpKa'])
    #os.system('eog %s' %file)
    #
    # Any difference for active site targets?
    #
    
    #
    # Plot it
    #
    nummuts={}
    nums=raw_data.keys()
    nums.sort()
    for num in nums:
        for co in raw_data[num].keys():
            max_val=-1.0
            sum=0.0
            count=0
            for dpka in raw_data[num][co]:
               if abs(dpka)>max_val:
                    max_val=abs(dpka)
               if dpka>0.01:
                   sum=sum+abs(dpka)
                   count=count+1
               #
               # Sort as function of number of mutations for other stat
               #
               if not nummuts.has_key(num):
                   nummuts[num]=[]
               nummuts[num].append(abs(dpka))
            if count==0:
                raw_data[num][co]=0
            else:
                raw_data[num][co]=float(sum)/float(count)
            #raw_data[num][co]=max_val
    import dislin_driver
    #dislin_driver.colour_2D(raw_data,'','','# of mutations','distance from target (A)','abs(dpKa)','dpka.tif')
    import os
    #os.system('gimp dpka.tif')
    #
    # Get dpKa as a function of # of mutants
    #
    #nums=nummuts.keys()
    #nums.sort()
    #x=[]
    #y=[]
    #for num in nums:
    #    for dpka in nummuts[num]:
    #        x.append(num)
    #        y.append(dpka)
    #file=dislin_driver.graf_mult2(x,[y],
    #                              title='dpKa, number of mutations',
    #                              x_legend='Number of mutations',
    #                              y_legend='abs(dpKa)')
    #os.system('gimp %s' %file)
    #
    # Save bigdict
    #
    fd=open('/home/nielsen/pKa-design/done_distnummuts/bigdict','w')
    import pickle
    pickle.dump(big_dict,fd)
    fd.close()

if __name__=='__main__':
    main()

