#!/usr/bin/env python

import os
PDBdir=os.path.join(os.getcwd(),'PDBs')
if not os.path.isdir(PDBdir):
    os.mkdir(PDBdir)


            

class pKa_contact_order:
    """Class for calculating the pKa contact order of proteins"""
    
    def __init__(self,pdbfiles):
        for pdbfile in pdbfiles:
            PI=self.get_PDB(pdbfile)
            RCOs=[PI.calculate_RCO(ss=False),PI.calculate_RCO(ss=True)]
            residues=sorted(RCOs[0].keys())
            residues.remove('all')
            print '%10s %7s %7s %7s' %('Residue','CO_noss','CO_w_ss','diff')
            for residue in residues:
                print '%10s %4s %7.2f %7.2f %7.2f' %(residue,PI.resname(residue),RCOs[0][residue],RCOs[1][residue],RCOs[1][residue]-RCOs[0][residue])
            print '------------------'
            print '%10s %4s %7.2f %7.2f %7.2f' %('All',' ',RCOs[0]['all'],RCOs[1]['all'],RCOs[1]['all']-RCOs[0]['all'])
            
        return
        
    #
    # -----
    #
            
    def get_PDB(self,name):
        """Get a PDB chain given as <PDBID><chainname>"""
        import Protool
        X=Protool.structureIO()
        CID=name[-1]
        name=name[:-1]
        #
        import os
        pdbfile=os.path.join(PDBdir,name+'.pdb')
        if not os.path.isfile(pdbfile):
            print 'Getting %s from the PDB' %(os.path.split(pdbfile)[-1])
            import Protool.PDBServices as PS
            PDB=PS.PDBServices()
            lines=PDB.getPDB(name)
            if len(lines)<100:
                import string
                raise Exception('PDB file not found: %s' %(string.join(lines)))
            X.parsepdb(lines)
            X.writepdb(pdbfile)
        else:
            X.readpdb(pdbfile)
        X.RemoveALT() # Removes alternative atoms
        X.Remove_All_NonAminoAcids() # Deletes all ligands and waters
        #
        # Keep only the chain we are interested in
        #
        for residue in X.residues:
            thisCID=X.chainid(residue)
            if CID=='_' and thisCID=='':
                pass
            elif CID!='_' and thisCID==CID:
                pass
            else:
                X.remove_residue(residue)
        X.Update()
        return X
        
class unfolded_pkas(pKa_contact_order):
    
    def __init__(self,pdbfiles):
        """Calculate a score for sequence-local contacts"""
        import Protool
        self.X=Protool.structure()
        
        
        for pdbfile in pdbfiles:
            xs=[]
            ys=[]
            labels=[]
            count=0
            PI=self.get_PDB(pdbfile)
            seq=PI.PirSeq(PI.sequence)
            window=4*[None]
            residues=4*[None]
            for residue,type in PI.sequence[:3]:
                window.append(type)
                residues.append(residue)
            for residue,type in PI.sequence[3:]:
                
                window.append(type)
                window=window[1:]
                #
                residues.append(residue)
                residues=residues[1:]
                #print window
                if window[3] in ['ASP','GLU']: #,'LYS','HIS','ARG']:
                    #print window
                    score=self.score_window(window)
                    print residues[3],window[3],score
                    xs.append(count)
                    labels.append(residues[3]+':%s' %window[3])
                    ys.append(score)
                    count=count+1
            #
            # Plot a bar chart
            #
            import pylab
            pylab.bar(xs,ys)
            import numpy
            xs=numpy.array(xs)
            pylab.xticks(xs+0.5,labels,rotation='vertical',size='x-small')
            pylab.title('dpKa values')
            pylab.show()
        return
        
    #
    # ----
    #
                
    def score_window(self,window):
        size=len(window)
        exclude=(size-1)/2
        center_res=window[exclude]
        if center_res in ['ASP','GLU']:
            gamma=-1
        elif center_res in ['HIS','LYS','ARG']:
            gamma=1
        else:
            return 0.0
        #
        score=0.0
        for count in range(size):
            dist=abs(exclude-count)
            res=window[count]
            if not res:
                continue
            if count==exclude:
                continue
            res=self.X.three_to_one[res]
            if res in 'DE':
                score=score+1.0/(dist**2)
            elif res in 'RKH':
                score=score-1.0/(dist**2)
            elif res in 'STYNQ':
                score=score+0.1/(dist**6)
            #elif res in 'FILWV':
            #    score=score-0.2*gamma/(dist**6)
        return score
            
        
if __name__=='__main__':
    #Y=unfolded_pkas(['2LZTA'])
    X=pKa_contact_order(['2LZTA'])
