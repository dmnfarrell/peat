#!/usr/bin/env python
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

#

class double_mutant:

    def doublemuts(self,stabdata,pKa_values):
        """Find all double mutant cycles and analyze them as a function of pH"""
        #
        # Standardize the names
        #

        for key in stabdata.keys():
            muts=key.split('+')
            import string
            muts=string.join(sorted(muts),'+')
            stabdata[muts]=stabdata[key]
        #
        dbl_cycles=self.find_double_mutants(stabdata.keys())
        for cycle in dbl_cycles:
            self.get_coupling(cycle,stabdata)
        return
        
    #
    # ------
    #
    
    def get_coupling(self,cycle,stabdata):
        """Get the coupling energy as a function of pH for this cycle"""
        print stabdata.keys()
        print 'Looking at cycle ',cycle
        count=0
        xs=[]
        ys=[]
        for ph in stabdata['ph']:
            Ewt=stabdata[cycle[0]][count]
            Emut1=stabdata[cycle[1]][count]
            Emut2=stabdata[cycle[2]][count]
            Edbl=stabdata[cycle[3]][count]
            cpl1=(Edbl-Emut1)-(Emut2-Ewt)
            cpl2=(Edbl-Emut2)-(Emut1-Ewt)
            print 'Coupling energy at pH %5.1f is %5.1f and %5.1f kJ/mol' %(float(ph),cpl1,cpl2)
            xs.append(float(ph))
            ys.append(cpl1)
            count=count+1
        import pylab
        pylab.plot(xs,ys,'ro-',label=str(cycle[4]))
        pylab.xlabel('pH')
        pylab.ylabel('interaction energy (kJ/mol)')
        pylab.legend()
        pylab.savefig('IEs/%s' %str(cycle[4]),dpi=300)
        pylab.clf()
        return
        
        
    #
    # ------
    #

    def addmuts(self,mutlist):
        """Add mutations"""
        import string
        muttxt=string.join(mutlist,'+')
        wt=False
        mutdef=[]
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
                        if not mut in sorted_muts:
                            allfound=False
                    #
                    # Everything is ok. Now which interaction are we measuring?
                    #
                    if allfound:
                        dbl_muts.append([wt.replace('+wt',''),leg1.replace('+wt',''),leg2.replace('+wt',''),dbl.replace('+wt',''),
                                        sorted([mut2.replace('+wt',''),mut3.replace('+wt',''),'Context: %s' %mut1])])
        return dbl_muts
        
        
        