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

import sys
import copy
import os

class wrapper_analyser:
    def go(self,file):
        try:
            lines=open(file).readlines()
            self.file = file
            os.system('mkdir %s' %self.file[:-4])
        except:
            print 'Error opening file',file
            return

        self.contents = {1:'Error tit',2:'Error activity',3:'Total error'}
            
        self.configs= {}

        for i in range(1,len(lines)):
            # find number of residues
            words=lines[i].split()
            incl_res=[]

            m = 4
            for w in words:
                try:
                    int(w)
                    break
                except:
                    m+=1
                
            for w in words[:m]:
                incl_res.append(w)
        
            t_incl_res = tuple(incl_res)

            if not self.configs.has_key(t_incl_res):
                self.configs[t_incl_res]={}

            # insert results
            c=1
            for w in words[len(t_incl_res):]:
                self.configs[t_incl_res][self.contents[c]] = w
                c+=1
            
                 
        res = self.configs.keys()
        res.sort(lambda x,y:self.sort(x,y))

        res_no_activity = []
        res_activity = []

        for r in res:
            if r[-1]=='-':
                res_no_activity.append(r)
            else:
                res_activity.append(r)
                
                #        print 'res_no_activit',res_no_activity
                #        print 'res_activit',res_activity
        self.comp_fits_no_act(res_no_activity)
        self.comp_fits_act(res_activity)

        return 
    


    def comp_fits_no_act(self,res):   
        # Compare fits without pH-activity profile
        curres = ''
        compare_error=0
        best_improvement = 0
        best_setup = ''
        out = ''
        all = open(self.file[:-4]+'/'+self.file[:-4]+'_no_activity.txt','w')
        all.write('Loaded residues      Best setup\n')
        for r in res:
            if not curres == self.residues(r):
                if not best_setup=='':
                    all.write('%-20s %s' %(str(curres),str(best_setup))+'\n')
                # new setup
                best_improvement = 0
                curres = self.residues(r)
                compare_error = float(self.configs[r][self.contents[1]])
                print 'Now comparing to',curres
                try:
                    out.close()
                except:
                    pass
                fn = self.file[:-4]+'/'+self.file[:-4]+'_no_activity_'+curres+'.txt'
                out = open(fn,'w')
                header = 'Residues, no. groups, no. acids, no. bases, no. deprotonated, Error on tits, Improvement\n'
                out.write(header)
                suffix = 'C\n'
                this_improvement = 0.0
                best_setup = r
            else:
                # compare
                this_improvement = (compare_error -
                                    float(self.configs[r][self.contents[1]]))/compare_error
                if this_improvement - best_improvement > 0.05:
                    best_improvement = this_improvement
                    best_setup = r
                    suffix = '*\n'
                else:
                    suffix = '\n'

            u=''
            for i in r:
                u+=i+' '
            u+=' %s'%self.configs[r][self.contents[1]]
            u+=' %6.2f%%'%(this_improvement*100)
            out.write(u+suffix)
        all.write('%-20s %s' %(str(curres),str(best_setup))+'\n')
        all.close()
        return


    def comp_fits_act(self,res):   
        # Compare fits that include pH-activity profile
        curres = ''
        compare_error_tit=0
        compare_error_act=0
        compare_error_tot=0
        best_improvement_tit = 0
        best_improvement_act = 0
        best_improvement_tot = 0
        out = ''
        best_setup_tit = ''
        best_setup_act = ''
        best_setup_tot = ''
        
        all = open(self.file[:-4]+'/'+self.file[:-4]+'_activity.txt','w')
        all.write('Loaded residues      Best setup for tits   Best setup for activity   Overall best setup\n')
        
        for r in res:
            if not curres == self.residues(r):
                if not best_setup_act=='':
                    all.write('%-20s %s %s %s' %(str(curres),str(best_setup_tit),str(best_setup_act),str(best_setup_tot))+'\n')
                # new setup
                best_improvement_tit = 0
                best_improvement_act = 0
                best_improvement_tot = 0

                curres = self.residues(r)
                compare_error_tit = float(self.configs[r][self.contents[1]])
                compare_error_act = float(self.configs[r][self.contents[2]])
                compare_error_tot = (float(self.configs[r][self.contents[1]])+float(self.configs[r][self.contents[2]]))/2.0
                print 'Now comparing to (incl activity data) ',curres
                try:
                    out.close()
                except:
                    pass
                fn = self.file[:-4]+'/'+self.file[:-4]+'_activity_'+curres+'.txt'
                out = open(fn,'w')
                header = 'Residues, no. groups, no. acids, no. bases, no. deprotonated, Error on tits, Error on activity data, Total error, Improvement on tits, Improvement on act, Total improvement\n'
                out.write(header)
                suffix = 'CCC\n'
                this_improvement_tit = 0.0
                this_improvement_act = 0.0
                this_improvement_tot = 0.0
                best_setup_tit = r
                best_setup_act = r
                best_setup_tot = r
                
            else:
                # compare
                this_improvement_tit = (compare_error_tit -
                                        float(self.configs[r][self.contents[1]]))/compare_error_tit
                this_improvement_act = (compare_error_act -
                                        float(self.configs[r][self.contents[2]]))/compare_error_act
                this_improvement_tot = (compare_error_tot -
                                        (float(self.configs[r][self.contents[1]])+
                                         float(self.configs[r][self.contents[2]]))/2.0
                                         )/compare_error_tot

                suffix = ''

                if (this_improvement_tit - best_improvement_tit) > 0.05:
                    best_improvement_tit = this_improvement_tit
                    best_setup_tit = r
                    suffix += '*'
                else:
                    suffix += ' '
                    
                if (this_improvement_act - best_improvement_act) > 0.05:
                    best_improvement_act = this_improvement_act
                    best_setup_act = r
                    suffix += '*'
                else:
                    suffix+=' '
                
                if (this_improvement_tot - best_improvement_tot) > 0.05:
                    best_improvement_tot = this_improvement_tot
                    best_setup_tot = r
                    suffix += '*\n'
                else:
                    suffix+=' \n'
                    

            u=''
            for i in r:
                u+=i+' '
            u+=' %7f'%float(self.configs[r][self.contents[1]])
            u+=' %7f'%float(self.configs[r][self.contents[2]])
            u+=' %7f'%((float(self.configs[r][self.contents[1]])+
                     float(self.configs[r][self.contents[2]]))/2.0)
            u+=' %6.2f%%'%(this_improvement_tit*100)
            u+=' %6.2f%%'%(this_improvement_act*100)
            u+=' %6.2f%%'%(this_improvement_tot*100)
            out.write(u+suffix)
        all.write('%-20s %s %s %s' %(str(curres),str(best_setup_tit),str(best_setup_act),str(best_setup_tot))+'\n')
        all.close()

        return






    def sort(self, x, y):
        # sort by no residues
        if not len(x) == len(y):
            return cmp(len(x),len(y))
        # sort by residues
        elif not self.residues(x) == self.residues(y):
            return cmp(self.residues(x), self.residues(y))
        # sort by no groups
        elif not x[-4] == y[-4]:
            return cmp(x[-4], y[-4])
        # sort by no acids
        elif not x[-3] == y[-3]:
            return cmp(y[-3], x[-3])
        # sort by no deprotonated
        else:
            return cmp(x[-1], y[-1])

    def residues(self,x):
        res = ''
        for w in x[:-4]:
            res+=w[:-6]
        return res


if __name__ == '__main__':
    wa = wrapper_analyser()
    wa.go(sys.argv[1])

