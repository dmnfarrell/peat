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

def get_net_charge(pdbfile,HIS):
    """Get the net charge within 20 A of the HIS"""
    import Protool
    X=Protool.structureIO()
    X.readpdb(pdbfile)
    close=[]
    HIS_ND1='%s:ND1' %HIS
    HIS_NE2='%s:NE2' %HIS
    for residue in X.residues.keys():
        for atom in X.residues[residue]:
            #print atom
            mdist=min(X.dist(HIS_ND1,atom),X.dist(HIS_NE2,atom))
            if mdist<50.0:
                close.append(residue)
                break
            elif mdist>355.0:
                break
    # Got all close residues, now count charge
    charge=0.0
    nc={'ASP':-1,'GLU':-1,'LYS':+1,'ARG':+1,'HIS':+1}
    close.sort()
    print close
    for res in close:
        restype=X.resname(res)
        if nc.has_key(restype):
            charge=charge+nc[restype]
            print res,restype,nc[restype],charge
    print 'Net charge',charge
    return charge

def get_HIS(pdbfilename,PD,Nuc):
    #
    # find the residue number of the HIS
    #
    import Protool
    X=Protool.structureIO()
    X.readpdb(pdbfilename)
    dists=[]
    for residue in X.residues.keys():
        if X.resname(residue)=='HIS':
            for PD_atom in ['OE1','OE2']:
                for HIS_atom in ['ND1','NE2']:
                    dist=X.distance(PD+':%s' %PD_atom,residue+':%s' %HIS_atom)
                    dists.append([dist,residue,PD+':%s' %PD_atom,residue+':%s' %HIS_atom])
    dists.sort()
    if dists[0][0]<3.5:
        print 'HIS Hbonded to PD: %s  %5.2f' %(dists[0][1],dists[0][0])
        HIS=dists[0][1]
        HIS_atom=dists[0][3][-3:]
        #print dists
        #print HIS_atom
        if HIS_atom=='ND1':
            other_atom=dists[0][3][:-3]+'NE2'
        elif HIS_atom=='NE2':
            other_atom=dists[0][3][:-3]+'ND1'
        else:
            raise Exception('Something is rotten')
        print 'Other atom',other_atom
        #
        # find the other Hbond partner of this HIS
        #
        HisDists=[]
        for residue in X.residues.keys():
            if residue==HIS or residue==PD:
                continue
            for atom in X.residues[residue]:
                dist=X.distance(other_atom,atom)
                if dist>30.0:
                    break
                HisDists.append([dist,atom])
        HisDists.sort()
        Hbond_res=[]
        for dist,atom in HisDists:
            if dist<3.2:
                Hbond_res.append([dist,atom])
        if len(Hbond_res)==0:
            Hbond_res.append(HisDists[0])
        print Hbond_res
        if len(Hbond_res)>2:
            raise Exception('More than two Hbond donors')
        other_side=[]
        for dist,atom in Hbond_res:
            print 'The residue "on the other side" is %s %s. dist: %5.2f A' %(atom,X.resname(atom),dist)
            other_side.append([pdbfilename,X.resname(atom)])
        #
        # Can we say anything about the protonation state of the HIS?
        #
        pass
        

    elif pdbfilename.find('H241A'):
        return 'Mutated',[]
    else:
        print 'No HIS Hbond found'
        print dists
        HIS=None
        stop
    return HIS,other_side
    

def main():
    import os
    scores={}
    for calctype in ['APO','HOLO','Helens_calcs']:
        dirs=os.listdir(calctype)
        for direc in dirs:
            if direc in ['CCPS','CCPS_3','CCPS_HIS+','CCPS_HIS-']:
                continue
            rdir=os.path.join(os.getcwd(),calctype,direc)
            if not os.path.isdir(rdir):
                continue
            files=os.listdir(rdir)
            found=False
            for fn in files:
                if fn.find('sensresult_10')!=-1:
                    found=fn
            
            print os.path.join(calctype,direc)
            name=direc
            if scores.has_key(direc):
                raise Exception('Duplicate PDBID!!: %s' %direc)
            #
            if found:
                fd=open(os.path.join(rdir,found))
                import pickle
                data=pickle.load(fd)
                fd.close()
                pdbfile=data[0]
                SR=data[1]
                PD=data[2]
                Nuc=data[3]
                #
                # get the Histidine pKa, and figure out if the His should be positive in the CCPS
                #
                if PD and Nuc:
                    HIS,other_side=get_HIS(os.path.join(rdir,pdbfile),PD[:-4],Nuc[:-4])
                    #print 'HIS',HIS
                #
                # Get the net charge within 20 A
                #
                if HIS.lower()!='mutated' and HIS:
                    net_charge=get_net_charge(os.path.join(rdir,pdbfile),HIS)
                else:
                    net_charge=None
                #
                # Do sensitivity analysis
                #
                import pKaTool.Do_Sensitivity_analysis as dosens
                score,csvline=dosens.print_result(pdbfile,SR,PD,Nuc,resultdir=rdir)
                scores[direc]={'score':score,'HIS':HIS,'HisPartner':other_side,'net_charge':net_charge}
                
                scores[direc]['PD_pKa']=dosens.get_pKa(PD,resultdir=rdir,calc=pdbfile)
                scores[direc]['Nuc_pKa']=dosens.get_pKa(Nuc,resultdir=rdir,calc=pdbfile)
                scores[direc]['His_pKa']=dosens.get_pKa(HIS,resultdir=rdir,calc=pdbfile)
                #
                # Get intrinsic pKa values
                #
                scores[direc]['PD_intpKa']=dosens.get_intpKa(PD,resultdir=rdir,calc=pdbfile)
                scores[direc]['Nuc_intpKa']=dosens.get_intpKa(Nuc,resultdir=rdir,calc=pdbfile)
                scores[direc]['His_intpKa']=dosens.get_intpKa(HIS,resultdir=rdir,calc=pdbfile)
                
            else:
                scores[direc]={'score':'notdone'}
                print 'No sensitivity analysis found'
                print '---------------------------------'
            if direc.find('ligand')!=-1 or direc.find('lignad')!=-1:
                ligand=True
            else:
                ligand=False
            scores[direc]['ligand']=ligand
            scores[direc]['csvline']=csvline
    #
    # Save the scores
    #
    fd=open('scores.pickle','w')
    import pickle
    pickle.dump(scores,fd)
    fd.close()
    return
    
    
def make_plots():
    fd=open('scores.pickle')
    import pickle
    scores=pickle.load(fd)
    fd.close()
    #
    # Do final stats
    #
    ligpdbs=[]
    for pdb in scores.keys():
        if scores[pdb]['ligand']:
            ligpdbs.append(pdb[:4])

    print 'Final stats'

    for structure in ['true_apo','apo_generated','holo']:
        notdone=0
        missingpka=0
        correct=0
        corrects=[]
        wrong=0
        wrongs=[]
        inconclusive=0
        incon=[]
        totcalcs=0
        printlines=[]
        partners={'normal':[],'reversed':[],'incon':[]}
        print 'Type of structures',structure
        PD_pKas=[]
        Nuc_pKas=[]
        His_pKas=[]
        PD_intpKas=[]
        Nuc_intpKas=[]
        His_intpKas=[]
        nc=[]
        cor_nc=[]
        rev_nc=[]
        inc_nc=[]
        for pdb in scores.keys():
            score=scores[pdb]['score']
            partner=None
            if scores[pdb].has_key('HisPartner'):
                partner=scores[pdb]['HisPartner']
            # Look only at the correct set of structures
            if structure=='true_apo':
                if pdb[:4] in ligpdbs or scores[pdb]['ligand']:
                    continue
            elif structure=='apo_generated':
                if not pdb[:4] in ligpdbs or scores[pdb]['ligand']:
                    continue
            elif structure=='holo':
                if not scores[pdb]['ligand']:
                    continue
                    
            if scores[pdb].has_key('csvline'):
                printlines.append(scores[pdb]['csvline'])
            if scores[pdb].has_key('net_charge'):
                if scores[pdb]['net_charge']:
                    nc.append(scores[pdb]['net_charge'])
            #
            # Keep track of the pKas
            #
            if score!='notdone':
                #print scores[pdb]
                for name,l1,l2 in [['PD_pKa',PD_pKas,PD_intpKas],
                                ['Nuc_pKa',Nuc_pKas,Nuc_intpKas],
                                ['His_pKa',His_pKas,His_intpKas]]:
                                pKa=scores[pdb][name]
                                intpKa=scores[pdb][name.replace('pKa','intpKa')]
                                if not pKa is None:
                                    l1.append(pKa)
                                    l2.append(intpKa)
                #PD_pKas.append(scores[pdb]['PD_pKa'])
                #Nuc_pKas.append(scores[pdb]['Nuc_pKa'])
                #His_pKas.append(scores[pdb]['His_pKa'])
            #
            # Tally the scores
            #
            totcalcs=totcalcs+1
            if score is None:
                missingpka=missingpka+1
                print 'MISSING',pdb
                raw_input('kkk')
            elif score == 'notdone':
                notdone=notdone+1
            elif score>0.0:
                correct=correct+1
                corrects.append(pdb)
                partners['normal'].append(partner)
                if scores[pdb].has_key('net_charge'):
                    if scores[pdb]['net_charge']:
                        cor_nc.append(scores[pdb]['net_charge'])
            elif score<0.0:
                wrong=wrong+1
                wrongs.append(pdb)
                partners['reversed'].append(partner)
                if scores[pdb].has_key('net_charge'):
                    if scores[pdb]['net_charge']:
                        rev_nc.append(scores[pdb]['net_charge'])
            elif score==0:
                inconclusive=inconclusive+1
                partners['incon'].append(partner)
                if scores[pdb].has_key('net_charge'):
                    if scores[pdb]['net_charge']:
                        inc_nc.append(scores[pdb]['net_charge'])
                # If inconclusive, then what is the other residue?
        printlines.sort()
        for line in printlines:
            print line

        wrongs.sort()
        corrects.sort()
        incon.sort()
        print 'Total calcs : %d' %totcalcs
        print 'Calcs not analyzed: %d' %notdone
        print 'Calcs analyzed: %d' %(totcalcs-notdone)
        print 'Correct IDs : %d' %correct
        print 'Reverse IDs   : %d' %wrong
        print 'Inconclusive: %d' %inconclusive
        print 'Missing pKas: %d' %missingpka
        #print 'Correct pdbs',corrects
        #print 'Wrong pdbs',wrongs
        print 'inconclusive'
        for ptype in partners.keys():
            print ptype
            print partners[ptype]
            print
        print '---------------------------'
        print
        import pylab
        
        print cor_nc
        print rev_nc
        print inc_nc
        n,bins,patches=pylab.hist([cor_nc,rev_nc,inc_nc],histtype='bar')
        
        print n
        print bins
        print patches
        #pylab.hist(rev_nc,label='reverse',histtype='barstacked')
        #pylab.hist(inc_nc,label='inconclusive',histtype='barstacked')
        pylab.legend([patches[0][0],patches[1][0],patches[2][0]],['normal','reverse','outlier'])
        pylab.title(structure)
        #pylab.show()
        # Histograms of pKa values
        #pylab.hist(PD_pKas,20,label='PD')
        #pylab.hist(Nuc_pKas,20,label='Nuc')
        #pylab.hist(His_pKas,20,label='His')
        #pylab.title(structure)
        #pylab.legend()
        #pylab.show()
         # Effect of delec
        delec_PD=[]
        delec_Nuc=[]
        ID_PD=[]
        ID_PDint=[]
        ID_Nuc=[]
        ID_Nucint=[]
        for count in range(len(PD_pKas)):
            delec_PD.append(PD_pKas[count]-PD_intpKas[count])
            if PD_pKas[count]>5.0 and PD_pKas[count]-Nuc_pKas[count]>=1.5:
                ID_PD.append(PD_pKas[count])
                ID_PDint.append(PD_intpKas[count])
                
            elif Nuc_pKas[count]>5.0 and Nuc_pKas[count]-PD_pKas[count]>=1.5:
                ID_Nuc.append(Nuc_pKas[count])
                ID_Nucint.append(Nuc_intpKas[count])
        for count in range(len(Nuc_pKas)):
            delec_Nuc.append(Nuc_pKas[count]-Nuc_intpKas[count])
        #Plot of pKa value against intrinsic pKa value
        #pylab.scatter(PD_intpKas,PD_pKas,len(PD_pKas),'b',label='PD',alpha=0.4)
        #pylab.scatter(ID_PDint,ID_PD,len(ID_PDint),'b',label='PD ID',alpha=1)
        
        #pylab.scatter(Nuc_intpKas,Nuc_pKas,len(Nuc_pKas),'r',label='Nuc',alpha=0.4)
        #pylab.scatter(ID_Nucint,ID_Nuc,len(ID_Nuc),'r',label='Nuc',alpha=1.0)
        # Plot the ones with ID
        
        #pylab.ylabel('pKa value')
        #pylab.xlabel('Intrinsic pKa value')
        #pylab.plot([0,12],[0,12],c='k',alpha=0.2)
        #pylab.plot([0,12],[5,5],c='g')
        #pylab.title('pKa values vs. Intrinsic pKa values')
        #pylab.legend()
        #pylab.xlim(0.0,12.0)
        #pylab.ylim(0.0,12.0)

        #pylab.show()
       
        #pylab.scatter(delec_PD,PD_pKas,len(PD_pKas),'b',label='PD',alpha=0.4)
        #pylab.scatter(delec_Nuc,Nuc_pKas,len(Nuc_pKas),'r',label='Nuc',alpha=0.4)
        #pylab.ylabel('pKa value')
        #pylab.xlabel('dpKa(elec)')
        #pylab.plot([0,12],[0,12],c='r')
        #pylab.plot([0,12],[5,5],c='g')
        #pylab.title('pKa values vs. delec values')
        #pylab.legend()
        #pylab.show()
        
    
if __name__=='__main__':    
    #import pylab
    #pkas=[4,4,5,5,4,3,8,8,9,7]
    #pylab.scatter(pkas,pkas,len(pkas))
    #pylab.hist(pkas,20)
    #pylab.show()
    #main()
    make_plots()

