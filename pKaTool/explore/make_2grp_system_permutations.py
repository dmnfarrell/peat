#!/usr/bin/env python

class permutations:

    def __init__(self):
        """Make permutations of a two-group acid-acid system"""
        self.pHstart=0.0
        self.pHend=20.0
        self.pHstep=0.1
        self.pKa_MCsteps=2000000
        
        import numpy
        pka_diffs=numpy.arange(0.0,2.0,0.01)
        int_enes=numpy.arange(0.0,2.0,0.01)
        #
        import math
        matrix={}
        for pka_diff in pka_diffs:
            print pka_diff
            matrix[pka_diff]={}
            for int_ene in int_enes:
                pKa_values, prot_states=self.do_pKa_calc(pka_diff,int_ene)
                #print int_ene,pKa_values
                #del prot_states[':0002:GLU']
                matrix[pka_diff][int_ene]=self.fit_curves(prot_states)
            #stop
        import TwoDplots
        TwoDplots.heatmap(matrix,'Two-group scan',firstkey='Difference in intrinsic pKa',secondkey='Interaction energy (kT)',zlabel='Average slope', 
        firstticks=[numpy.arange(0,200,10),numpy.arange(0.0,2.0,0.1)], 
        secondticks=[numpy.arange(0,200,10),numpy.arange(0.0,2.0,0.1)])
        #
        return
        
    #
    # -----
    #
    
    def fit_curves(self,curves):
        error=[]
        import string, pKarun
        PKana=pKarun.pKa_general.pKanalyse()
        for group in curves.keys():
            solution,sq=PKana.fit_to_henderson(curves[group])
            sq=solution[0]
            try:
                sq=float(sq)
            except:
                sq=0.0
            error.append(abs(sq))
        return sum(error)/len(error)
        
    #
    # -----
    #
    
    def do_pKa_calc(self,pka_diff,int_ene): 
        """Do the pKa calculation for this system"""
        import pKaTool.pKa_calc, math
        self.PK=pKaTool.pKa_calc.Boltzmann()
        groups={':0001:GLU':{'backgr':0.0,':0002:GLU':int_ene},':0002:GLU':{'backgr':pka_diff*math.log(10),':0001:GLU':int_ene}}
        self.setup_system(groups)
        #
        pKa_values,prot_states=self.PK.calc_pKas(mcsteps=self.pKa_MCsteps,
                                                 phstep=self.pHstep,
                                                 phstart=self.pHstart,
                                                 phend=self.pHend,
                                                 exp_pHs=[],
                                                 verbose=1,
                                                 complete_pka=False)
        return pKa_values,prot_states
        
        
        
    def setup_system(self,groups):
        """Setup the MC_CPP system for the nontit system"""
        # 
        # Get the intrinsic pKas
        #
        intpkas={}
        backgr={}
        desolv={}
        import math
        for group in groups.keys():
                intpkas[group]=0.0#groups[group]['intpka']
                backgr[group]=0.0
                desolv[group]=groups[group]['backgr']
        self.PK.groups=intpkas.keys()
        self.PK.groups.sort()
        self.PK.desolv=desolv.copy()
        self.PK.backgr=backgr.copy()
        #
        matrix={}
        prune_matrix={}
        for group1 in groups.keys():
            prune_matrix[group1]={}
            for group2 in groups.keys():
                if groups[group1].has_key(group2):
                    prune_matrix[group1][group2]=[groups[group1][group2],0.0,0.0,0.0]
                else:
                    prune_matrix[group1][group2]=[0.0,0.0,0.0,0.0]
        self.PK.matrix=prune_matrix.copy()
        return
        
        
if __name__=='__main__':
    X=permutations()