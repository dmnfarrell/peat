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
ln10=math.log(10)

def getfile(file):
    import cPickle
    fd=open(file)
    d=cPickle.load(fd)
    fd.close()
    return d

#
# -----
#

def get_files(only=None,type=None):
    #
    # Get the names of all the files
    #
    loadfiles=[]
    import os
    dir='/home/nielsen/pKa-design/done_distnummuts'
    dirs=os.listdir(dir)
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
            for filename in files:
                if filename[-7:]=='pka.pdb' or filename.find('8tln.pdb')!=-1:
                    loadfiles.append(os.path.join(rdir,filename))
    return loadfiles
#
# ----
#

def get_phidpka(target=None,mutation=None,su=None,mc=None,matrix=None):
    """Get the delta pKa from phi/ln(10) and the real dpKa from the MC calcs"""
    if su.has_key(mutation):
        #
        # This is the simple case. We have an interaction energy for the mutant,
        # so we can calculate the dpka easily
        #
        phidpka=su[mutation][target]/ln10
    else:
        #
        # If the mutant is not in the sugelm data then it could
        # be a neutral mutation
        #
        s=mutation.split(':')
        orgname=':'+s[1]+':'+s[0]
        phidpka=-matrix[target][orgname][0]/ln10
    #
    # If we have a base as the target then we expect the opposite pKa shift
    #
    import pKarun
    if not pKarun.isacid(target):
        phidpka=-phidpka
    #
    # Get the real observed shift
    #
    rdpka=mc[mutation][target]

    return phidpka,rdpka

#
# ----
#

def get_effective_eps(target,mutation,dpKa,X_mut=None,X_wt=None,phidpka=None):
    """Get the effective dielectric constant for a mutation"""
    if target.find('NTERM')!=-1 or target.find('CTERM')!=-1:
        return None,0
    if abs(dpKa)<0.1:
        return None,0
    #
    crg_group={'ASP':'CG',
               'GLU':'CD',
               'LYS':'NZ',
               'ARG':'CZ',
               'HIS':'CD2',
               'TYR':'OH',
               'GLN':'CD',
               'ASN':'CG',
               'PHE':'CG',
               'CYS':'SG',
               'PRO':'CG'}

    #
    # Is there are mutant PDB file?
    #
    dist=None
    if X_mut:
        #
        # Yes
        #
        target_num=X_mut.resnum(target)
        if X_mut.atoms.has_key(target_num+':'+crg_group[X_mut.resname(target)]):
            target_atom=target_num+':'+crg_group[X_mut.resname(target)]
            if not X_mut.atoms.has_key(target_atom):
                print 'could not find',target_atom
                print X_mut.residues[target_num]
            #
            # Find the mutated residue
            #
            mutated_num=':'+mutation.split(':')[1]
            mutated_resname=X_mut.resname(mutated_num)
            new_residue=mutation.split(':')[-1]
            if new_residue!=mutated_resname:
                print 'Mutated PDB file is corrupt!'
                print mutation
                raise 'corrupt mutated PDB file'
            mutated_atom=mutated_num+':'+crg_group[mutated_resname]
            dist=X_mut.distance(mutated_atom,target_atom)
        else:
            print 'did not find'
    else:
        #
        # No, - calculate distance from wt pdb file
        #
        target_num=X_wt.resnum(target)
        if X_wt.atoms.has_key(target_num+':'+crg_group[X_wt.resname(target)]):
            target_atom=target_num+':'+crg_group[X_wt.resname(target)]
            if not X_wt.atoms.has_key(target_atom):
                print 'could not find',target_atom
                print X_wt.residues[target_num]
            #
            # Find the mutated residue
            #
            mutated_num=':'+mutation.split(':')[1]
            mutated_resname=X_wt.resname(mutated_num)
            mutated_atom=mutated_num+':'+crg_group[mutated_resname]
            dist=X_wt.distance(mutated_atom,target_atom)
    #
    # Calculate the effective eps
    #
    if not dist:
        raise 'could not find mutation'
    #if abs(dpKa)>=0.2 and abs(phidpka-dpKa)<0.1:
    if abs(dpKa)>0.0001:
        ef_eps=243.40/(abs(dpKa)*dist)
    else:
        ef_eps=None
    if ef_eps:
        print '%15s %15s %5.2fA dpKa:%5.2f eps=%5.2f' %(target,mutation,dist,dpKa,ef_eps)
    else:
        print '%15s %15s %5.2fA dpKa:%5.2f eps=N/A' %(target,mutation,dist,dpKa)  

    return ef_eps,dist

#
# ----
#

def analyse_one_pdbfile(pdbfile,bigdict=None):
    """Load the MC.tabdata file and the sugelm file to determine the
    effective dielectric constant for single mutations"""
    mcfile=pdbfile+'.MC.tabdata'
    sugelm=pdbfile+'.sugelm_data'
    print 'Loading:'
    print mcfile
    print sugelm
    print
    import os
    if not os.path.isfile(mcfile):
        return [],[],0
    if not os.path.isfile(sugelm):
        return [],[],0
    if bigdict:
        fd=open(bigdict)
        import pickle
        big_dict=pickle.load(fd)
        fd.close()
    #
    # Get the PDB file
    #
    import pKaTool.pKaIO
    P=pKaTool.pKaIO.pKaIO(pdbfile)
    # Matrix
    matrix=P.read_matrix()
    # Intrinsic pKa values for mutated residues
    intpkafile=pdbfile+'.intpka_data'
    fd=open(intpkafile)
    import pickle
    mutant_intpkas=pickle.load(fd)
    fd.close()
    #
    # Read the PDB file
    #
    import Protool
    PDB=Protool.structureIO()
    PDB.readpdb(pdbfile)
    #
    # Start calculating
    #
    mc=getfile(mcfile)
    su=getfile(sugelm)['data']
    print 'Number of mutations in sugelm',len(su.keys())
    print 'Number of mutations in tabdata',len(mc.keys())

    sites={}
    import string
    for mutation in su.keys():
        orgres=mutation.split(':')[:2]
        orgres=string.join(orgres,':')
        sites[orgres]=1
    print 'Number of unique sites',len(sites.keys())
    #
    # Should we do the full solutions or just single mutations?
    #
    if bigdict:
        print 'Getting mutations from bigdict'
        mutations=[]
        for key in big_dict.keys():
            pdbfile_name=os.path.split(pdbfile)[1]
            pos=key.find(pdbfile_name)
            if pos!=-1:
                target=key[pos+len(pdbfile_name):]
                for muts,dpka in big_dict[key]:
                    mutations.append([target,muts,dpka])
    else:
        mutations=mc.keys()
    #
    # Load the wild type PDB file
    #
    X_wt=Protool.structureIO()
    X_wt.readpdb(pdbfile)
    #
    # go on, calculate the difference between the MC result and the phi result
    #
    mutations.sort()
    phi=[]
    real=[]
    ratio=[]
    dist=[]
    epses_dpKa=[]
    epses=[]
    import types
    for mutant in mutations:
        if type(mutant) is types.ListType:
            #
            # Deal with multiple mutations
            #
            sum_phidpka=0.0
            sum_MCdpka=0.0
            mutations=mutant[1]
            target=mutant[0]
            rdpka=mutant[2]
            if len(mutations)<2:
                continue
            for mutation in mutations:
                phidpka,MCdpka=get_phidpka(mutation=mutation,target=target,su=su,mc=mc,matrix=matrix)
                sum_phidpka=sum_phidpka+phidpka
                sum_MCdpka=sum_MCdpka+MCdpka
            phi.append(sum_phidpka)
            real.append(rdpka)
            #if abs(rdpka-sum_phidpka)>1.0:
            #    print 'T: %15s phidpka %5.1f rdpka: %5.1f muts: %s ' %(target,sum_phidpka,rdpka,str(mutations))
        else:
            #
            # This is for looking at single mutations
            #
            import Design_pKa_help
            X=Design_pKa_help.pKa_dist(pdbfile)
            targets=mc[mutant].keys()
            targets.sort()
            #
            # Load the mutant PDB file if we can
            #
            import os
            pdbdir=pdbfile+'.pdbs'
            mutfile=os.path.join(pdbdir,mutant+'.pdb')
            if os.path.isfile(mutfile):
                import Protool
                X_mut=Protool.structureIO_fast()
                X_mut.readpdb(mutfile)
            else:
                X_mut=None
            #
            # Loop over all targets
            #
            for target in targets:
                #
                # Get the distance between the target and the mutatn
                #
                targetnum=':'+target.split(':')[1]
                mutantnum=':'+mutant.split(':')[1]
                distance=X.get_min_dist(target,mutant)
                #
                # Get the delta pKa values
                #
                phidpka,rdpka=get_phidpka(mutation=mutant,target=target,su=su,mc=mc,matrix=matrix)
                if not rdpka:
                    continue
                if abs(phidpka)>=0.0001 and abs(rdpka)<20.0:
                    phi.append(phidpka)
                    real.append(abs(rdpka))
                    dist.append(distance)
                    #
                    # Effective eps
                    #
                    eps_eff,distance_eps=get_effective_eps(target,mutant,abs(rdpka),X_mut=X_mut,X_wt=X_wt,phidpka=phidpka)
                    if eps_eff:
                        epses.append(eps_eff)
                        #epses_dpKa.append(abs(rdpka))
                        epses_dpKa.append(distance_eps)
                    #ratio.append(rdpka/phidpka)
                    #print phidpka,rdpka
    tabdata_muts=len(mutations)
    return phi,real,tabdata_muts,dist,epses_dpKa,epses

#
# ----
#

if __name__=="__main__":
    files=get_files()
    xs_single=[]
    ys_single=[]
    xs_mult=[]
    ys_mult=[]
    dists=[]
    names=[]
    epses_all=[]
    epses_dpKa_all=[]
    muts=0
    for filename in files:
        print 'Processing',filename
        name=filename[-12:-7]
        names.append(name)
        x,y,nummuts,dist,epses_dpKa,epses=analyse_one_pdbfile(filename)
        xs_single.append(x)
        ys_single.append(y)
        dists.append(dist)
        epses_all.append(epses)
        epses_dpKa_all.append(epses_dpKa)
        print nummuts
        muts=muts+nummuts
        #
        # Now for multiple mutations
        #
        #x,y=analyse_one_pdbfile(filename,'/home/nielsen/pKa-design/done_distnummuts/bigdict')
        #xs_mult.append(x)
        #ys_mult.append(y)
    #
    # Plot real effective eps
    #
    import dislin_driver
    filename=dislin_driver.graf_mult3(xs=epses_dpKa_all,ys=epses_all,title='Effective dielectric constant for single mutations',legends=names,x_legend='distance (A)',y_legend='Effective dielectric constant')
    #
    # Plot it
    #
    count=0
    for x in xs_single:
        count=count+len(x)

    print 'Number of delta pKa values',count
    stop
    import dislin_driver
    #filename=dislin_driver.graf_mult3(xs=xs_single,ys=ys_single,title='pKa shifts for single mutations',legends=names,x_legend='dpKa (dphi/ln10)',y_legend='dpKa')
    #
    # Plot dpKa vs distance
    #
    filename=dislin_driver.graf_mult3(xs=dists,ys=ys_single,title='Distance dependence of dpKas',
                                       legends=names,
                                       x_legend='Distance (A)',
                                       y_legend='abs(dpKa)')
    #
    # Plot the multiple muts
    #
    #filename=dislin_driver.graf_mult3(xs=xs_mult,ys=ys_mult,
    #                                  title='Effective eps, multiple mutations',
    #                                  legends=names,
    #                                  x_legend='dpKa from phi/ln10',
    #                                  y_legend='dpKa from MC')
    print 'Total number of mutations',muts
    allx=[]
    for x in xs_single:
        allx=allx+x
    sum=0.0
    for x in allx:
        sum=sum+abs(x)
    print 'Average dpKa: %5.2f' %(sum/float(len(allx)))
