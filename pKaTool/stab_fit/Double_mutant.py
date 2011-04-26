#!/usr/bin/env python

class double_mutant:

    def doublemuts(self,stabdata,pKa_values):
        """Find all double mutant cycles and analyze them as a function of pH"""
        dbl_cycles=self.find_double_mutants(stabdata.keys())
        for cycle in dbl_cycles:
            print cycle
        stop
        return

    def addmuts(self,mutlist):
        """Add mutations"""
        import string
        muttxt=string.join(mutlist,'+')
        wt=False
        mutdef=[]
        print muttxt
        for mut in muttxt.split('+'):
            if mut=='wt':
                wt=True
            else:
                mutdef.append(mut)
        if wt:
            mutdef.append('wt')
        import string
        txt=sorted(mutdef)
        txt=string.join(txt,'+')
        print txt
        print '-=-=--=--=-'
        return txt
    
    #
    # -----
    #

    def find_double_mutants(self,mutants):
        """Find all double mutant cycles in a list of mutants"""
        dbl_muts=[]
        sorted_muts=[]
        for prot in mutants:
            import string
            muts=sorted(prot.split('+'))
            if not 'wt' in muts:
                muts.append('wt')
                muts.sort()
            sorted_muts.append(string.join(muts,'+'))
        #sorted_muts.append('')
        sorted_muts.remove('ph+wt')
        #
        # Now all mutants are sorted, then we try to split all possible ways
        #
        for mut1 in sorted_muts:
            for mut2 in sorted_muts:
                if mut2==mut1:
                    continue
                for mut3 in sorted_muts:
                    if mut3==mut1 or mut3==mut2:
                        continue
                    wt=mut1
                    leg1=self.addmuts([mut1,mut2])
                    leg2=self.addmuts([mut1,mut3])
                    dbl=self.addmuts([mut1,mut2,mut3])
                    # 
                    # Make sure all are different
                    #
                    cycle=[wt,leg1,leg2,dbl]
                    same=False
                    for pos1 in range(len(cycle)-1):
                        for pos2 in range(pos1+1,len(cycle)):
                            if cycle[pos1]==cycle[pos2]:
                                same=True
                    if same:
                        continue
                    #
                    # Check that all mutations exist
                    #
                    allfound=True
                    for mut in cycle:
                        print mut
                        if not mut in sorted_muts:
                            allfound=False
                    #
                    # Everything is ok. Now which interaction are we measuring?
                    #
                    print 'Measuring %s - %s' %(mut2,mut3)
                    if allfound:
                        dbl_muts.append([wt,leg1,leg2,dbl,sorted([mut2.replace('+wt',''),mut3.replace('+wt',''),'Context: %s' %mut1])])
                    print '============='
        return dbl_muts
        
        
        