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

from Tkinter import *
import pKa_system 
import os
import copy

class wrapper:
    def __init__(self,script=None):
        self.all_residues = []
        self.all_a_bs = []
        self.no_residues = 0
        self.residue_files = {}
        self.active_states = []
        self.runs=0
        self.dummy_run = 0       
        if script:
            lines = open(script,'r').readlines()
            self.interpret(lines)
        #else:

        #if script=='BCX_N35D':
        #    self.parameters_BCX_N35D()
        
        #if script=='BCX_WT':
        #    self.parameters_BCX_WT()

        #if script=='FN3':
        #    self.parameters_FN3()
        #if script=='ERT1':
        #    self.parameters_ERT1()
        #if script=='ERT2':
        #    self.parameters_ERT2()
        #if script=='ERT3':
        #    self.parameters_ERT3()
        #if script=='ERT4':
        #    self.parameters_ERT4()
        #self.go()


        return


##     def parameters_BCX_N35D(self):

##         self.open_report('/home/people/chresten/PhD/titration_data/BCX/id_no_groups_N35D/report_test_2.txt')
        
        
##         self.residue_files = {'D35':['/home/people/chresten/PhD/titration_data/BCX/bcx_n35d_asp_full.txt',4,8],
##                               'E78':['/home/people/chresten/PhD/titration_data/BCX/n35d_full.txt',7,7],
##                               'E172':['/home/people/chresten/PhD/titration_data/BCX/n35d_full.txt',6,7]}

##         self.ph_activity = '/home/people/chresten/PhD/titration_data/BCX/n35d_activity.txt'
##         self.include_ph = 1
##         return



##     def parameters_BCX_WT(self):

##         self.open_report('/home/people/chresten/PhD/titration_data/BCX/id_no_groups_WT/report_test_2.txt')
              
##         self.residue_files = {'E78':['/home/people/chresten/PhD/titration_data/BCX/wildtype_full.txt',7,7],
##                               'E172':['/home/people/chresten/PhD/titration_data/BCX/wildtype_full.txt',6,7]}

##         self.ph_activity = '/home/people/chresten/PhD/titration_data/BCX/wt_activity.txt'
##         self.include_ph = 1
##         return


##     def parameters_FN3(self):

##         self.open_report('/home/people/chresten/PhD/titration_data/FN3/id_no_groups/report_test_2.txt')
              
##         self.residue_files = {'Asp7':['/home/people/chresten/PhD/titration_data/FN3/Asp7.txt',1,1],
##                               'Glu9':['/home/people/chresten/PhD/titration_data/FN3/Glu9.txt',1,1],
##                               'Asp23':['/home/people/chresten/PhD/titration_data/FN3/Asp23.txt',1,1]}

##         self.include_ph = 0
##         return


##     def parameters_ERT1(self):

##         self.open_report('/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/id_no_group_dyson/report.txt')
              
##         self.residue_files = {'Asp26':['/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/Dyson1996/Asp26_lowpH.txt',1,1],
##                               'Cys32':['/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/Dyson1995/Cys32.txt',1,1],
##                               'Cys35':['/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/Dyson1995/Cys35.txt',1,1]}

##         self.ph_activity = '/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/activity_figure4.txt'
##         self.include_ph = 1
##         self.initial_pkas = 7.0
##         return

##     def parameters_ERT2(self):

##         self.open_report('/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/id_no_group_dyson/report2.txt')
              
##         self.residue_files = {'Asp26':['/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/Dyson1996/Asp26_all.txt',1,1],
##                               'Cys32':['/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/Chivers1997/Cys32.txt',1,1]}

##         self.ph_activity = '/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/activity_figure4.txt'
##         self.include_ph = 1
##         self.initial_pkas = 7.0
##         return

##     def parameters_ERT3(self):

##         self.open_report('/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/id_no_group_dyson/report3.txt')
              
##         self.residue_files = {'Asp26':['/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/Dyson1996/Asp26_all.txt',1,1],
##                               'Cys32':['/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/Chivers1997/Cys32.txt',1,1]}

##         self.ph_activity = '/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/activity_figure3.txt'
##         self.include_ph = 1
##         self.initial_pkas = 7.0
##         return

##     def parameters_ERT4(self):

##         self.open_report('/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/id_no_group_dyson/report4.txt')
              
##         self.residue_files = {'Asp26':['/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/Dyson1996/Asp26_lowpH.txt',1,1],
##                               'Cys32':['/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/Dyson1995/Cys32.txt',1,1],
##                               'Cys35':['/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/Dyson1995/Cys35.txt',1,1]}

##         self.ph_activity = '/home/people/chresten/PhD/titration_data/thioredoxins/EcoliReducedThioredoxin/activity_figure3.txt'
##         self.include_ph = 1
##         self.initial_pkas = 7.0
##         return



    def do_search(self):
        """Searches the loaded titration curves to find optimal fit
        with pH-activity profile"""

        # two titration curves (+ one unconstrained)
        for i in range(len(self.residue_files.keys())):
            for j in range(i+1,len(self.residue_files.keys())):
                r1=self.residue_files.keys()[i]
                r2=self.residue_files.keys()[j]

                print r1,r2
                self.run([r1,r2], ['a','a','1'])
                self.run([r1,r2], ['a','a','a','1'])
                    


        print '---------------------------------'
        print 'Done - Searched %d combinations' %self.runs


        return

        

    def all_combos(self):
        """Does fits to all combinations of acid/base assignment and number
        of deprotonated groups to less than five titration curves"""
        residues = self.residue_files.keys()

        # make combinations for one group
        self.generate_permutations(1,4,self.include_ph)

        for residue in residues:
            for p in self.permutations:
                #print p
                self.run([residue],p)
                
        # make combinations for two group
        self.generate_permutations(2,4,self.include_ph)

        for i in range(len(residues)):
            for j in range(i+1,len(residues)):
                for p in self.permutations:
                    #print p
                    self.run([residues[i],residues[j]],p)
        

        # make combinations for three group
        self.generate_permutations(3,4,self.include_ph)

        for p in self.permutations:
            #print p
            self.run(residues,p)
            
        self.close_report()
        print 'Total runs:',self.runs
        return


    def run(self, res, p):
        no_g = 0
        no_a = 0
        no_b = 0
        a_b = []
        for i in p:
            if i[0]=='a':
                a_b.append(1)
                no_g+=1
                no_a+=1
               
            if i[0]=='b':
                a_b.append(0)
                no_g+=1
                no_b+=1
   

        if p[-1]=='1':       
            states = self.generate_active_states(no_g,p[:-1])
        else:
            states=[[]]        

        for state in states:
            #
            # Launch for these residues with this permutation and this active state
            #
            print 'Do ',res,p,state
            self.launch(no_g)
            self.set_acid_base(a_b)
            
            for i in range(len(res)):
                self.residue = res[i]
                
                if p[i][0] == 'a':
                    self.a_b = 'acid'
                else:
                    self.a_b = 'base'
    
                self.load(self.residue_files[res[i]])

            if getattr(self,'initial_pkas',None):
                self.set_pkas(self.initial_pkas)

            if p[-1]=='1':
                self.load_ph_activity(self.ph_activity)
                self.set_CCPS(state[1])
                self.fit_to_curves_and_ph_activity_profile()
                self.report_errors(include_activity=1,no_deprotonated=state[0],no_groups=no_g,no_acids=no_a,no_bases=no_b)
                
            else:
                self.fit()
                self.report_errors(include_activity=0,no_groups=no_g,no_acids=no_a,no_bases=no_b)
        
            



            self.runs+=1
            self.reset()


    def generate_permutations(self, no_exp_groups, max_groups, include_ph):
        self.no_groups = no_exp_groups
        self.max_groups = max_groups
        self.frame = []

        for i in range(max_groups):
            self.frame.append(0)

        if include_ph:
            self.frame.append(0)


        self.permutations = []
        self.mutate(0,'a')

        return 


    def mutate(self,i,last):
        if i<self.no_groups:
            #combis =['a']
            combis =['aa','bb']
        elif(i<self.max_groups):
            if last in ['aa','bb','a']:
                combis = ['a','b','0']
            elif last=='b':
                combis = ['b','0']
            else:
                combis = ['0']
            
        elif(i<self.max_groups+1):
            combis=['0','1']
            
        for j in range(len(combis)):
            self.frame[i]=combis[j]
            if i<len(self.frame)-1:
                self.mutate(i+1,self.frame[i])
            else:
                temp = copy.copy(self.frame)
                self.permutations.append(temp)
        return 


    def generate_active_states(self,no_groups,p):

 #       print p


        #
        # build list of active states as they are listed in pKaTool
        #
        self.charged_states = []
        temp = []
        for n in range(no_groups):
            temp.append(0)
        self.make_charged_states(0,temp)


        states = []
        for no_deprotonated in range(no_groups):
 #           print 'looking for',no_deprotonated,'deprotonated'
            for j in range(len(self.charged_states)):
                cs = self.charged_states[j]
                this_no_deprot = 0
                for i in range(len(cs)):
                    if cs[i] == 1 and p[i][0]=='a':
                        this_no_deprot+=1
                    elif cs[i] == 0 and p[i][0]=='b':
                        this_no_deprot+=1
 #               print 'state',cs,'has',this_no_deprot,'deprotonated'
                if this_no_deprot == no_deprotonated:
                    temp2 = []
                    for k in range(len(self.charged_states)):
                        if not k==j:
                            temp2.append(0)
                        else:
                            temp2.append(1)
#                    print 'found',temp2
                    states.append([no_deprotonated,copy.deepcopy(temp2)])
                    break

        
#        states = {1:[[0,[1,0]],
#                     [1,[0,1]]],
#                  2:[[0,[1,0,0,0]],
#                     [1,[0,1,0,0]],
#                     [2,[0,0,0,1]]],
#                  3:[[0,[1,0,0,0,0,0,0,0]],
#                     [1,[0,1,0,0,0,0,0,0]],
#                     [2,[0,0,0,1,0,0,0,0]],
#                     [3,[0,0,0,0,0,0,0,1]]],
#                  4:[[0,[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],
#                     [1,[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],
#                     [2,[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0]],
#                     [3,[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]],
#                     [4,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]]]}
        

#        print states
        for s in states:
            print s
        return states
    



    def make_charged_states(self,no,temp):
        if not no == len(temp):
            for p in [0,1]:
                temp[no]=p
                self.make_charged_states(no+1,temp)
        else:
            self.charged_states.append(copy.deepcopy(temp))
           # print temp
        return
                




    def interpret(self,lines):
        for line in lines:
            if len(line) > 2:
                if not line[0:1] == '#':
                    print 'Got:',line[:-1]
                    self.command(line.split())
        return

    def command(self, words):
        if words[0]=='launch':
            self.launch(int(words[1]))
        elif words[0]=='report':        
            self.open_report(words[1])
        elif words[0]=='add_residue':
            self.add_residue_file(words[1], words[2], int(words[3]), int(words[4]))
        elif words[0]=='a/bs':
            self.a_b=words[1]
        elif words[0]=='loadtits':
            self.load(words[1])
        elif words[0]=='activegroups':
            self.set_active_groups(words[1:])
        elif words[0]=='fit':
            self.fit()
        elif words[0]=='loadactivity':
            self.load_ph_activity(words[1])
        elif words[0]=='set_activity_file':
            self.set_ph_activity(words[1])
        elif words[0]=='activestates':   
            self.set_CCPS(words[1:])
        elif words[0]=='fit_tits_ph':
            self.fit_to_curves_and_ph_activity_profile()
        elif words[0]=='reporterrors':
            self.report_errors(include_activity=int(words[1]))
        elif words[0]=='do_search':
            self.do_search()
        elif words[0]=='set_dummy_run':
            self.set_dummy_run(words[1])    
        elif words[0]=='all_combos':
            self.all_combos()


        return

    def set_dummy_run(self,value):
        self.dummy_run = int(value)
        return

    def add_residue_file(self, residue, file, pos, no_res_in_file):
        self.residue_files[residue] = [file,pos,no_res_in_file]
        return
    
    def launch(self,no_groups):
        self.pka = pKa_system.pKa_system(no_groups)
        return
    
    def open_report(self,file):
        self.report = open(file,'w')
        self.report.write('[Loaded groups], no_groups, no_deprotonated, error_tits, error_ph_actvity\n')

        return

    def close_report(self):
        self.report.close()
        return
    
    def load(self,filename):
        acid_bases = []
        self.all_residues.append(self.residue)
        self.all_a_bs.append(self.a_b)
        print filename
        for i in range(1,filename[2]+1):
            temp = ['', StringVar(),StringVar()]
            if i == filename[1]:
                temp[0]=self.residue
                temp[1].set(self.a_b)
                temp[2].set(':%04d:ASP'%(self.no_residues))
                print 'loading--------->',temp[2].get()
                self.no_residues+=1
            else:
                temp[2].set('None')

            acid_bases.append(temp)
            
        for ab in acid_bases:
            print ab[0],ab[1].get(),ab[2].get()

        self.pka.load_curves(filename=filename[0],acid_bases=acid_bases)
        return

    def set_pkas(self,pka):
        for group in self.pka.groups:
            self.pka.groups[group].intpka.set(pka)
            print 'setting pka ',pka

        self.pka.update_scales_from_fit()
        self.pka.titwin.update()
        return

    def reset(self):
        self.all_residues = []
        self.all_a_bs = []
        self.no_residues = 0
        self.active_states = []
        return

    def fit(self):
        if not self.dummy_run:
            self.pka.fit_system_to_curves()
        return

    def set_active_groups(self, actives):
        for i in range(len(actives)):
            self.pka.groups[i].active.set(int(actives[i]))
        return

    def set_acid_base(self, a_b):
        for i in range(len(a_b)):
            self.pka.groups[i].acid_base.set(int(a_b[i]))
        return

    def set_ph_activity(self,filename):
        print 'ph file',filename
        self.ph_activity = filename
        self.include_ph = 1
        return
    
    def load_ph_activity(self,filename):
        self.pka.load_pH_activity_profile(filename=filename)
        return

    def set_CCPS(self, actives):
        for i in range(len(actives)):
            self.pka.act_state[i].set(int(actives[i]))
            self.active_states.append(int(actives[i]))
        return
    
    def fit_to_curves_and_ph_activity_profile(self):
        if not self.dummy_run:
            self.pka.fit_to_curves_and_ph_activity_profile()
        return
    
    def get_errors(self,include_activity=0):
        errors = []
        # error of titration curves
        curve_error = self.pka.evaluate_fit(make_win=0)[1]
        errors.append(curve_error)
        #error of pH-activity profile
        if include_activity:
            activity_error = self.pka.evaluate_fit_ph_activity(make_win=0)
            errors.append(activity_error)

        return errors

    def report_errors(self,include_activity=0,no_deprotonated=None, no_groups=0, no_acids=None, no_bases=None):
        """Writes current errors to the report file"""

        # get errors of current fit
        errors = self.get_errors(include_activity=include_activity)

        # generate output line for report 
        line = ''
        # residues
        for i in range(len(self.all_residues)):
            line+= '%4s(%4s) '%(self.all_residues[i],self.all_a_bs[i])
        # number of groups
        line+='%2d '%no_groups
        # number of acids/bases
        if not no_acids==None:
            line+='%2d '%no_acids
        if not no_bases==None:
            line+='%2d '%no_bases
        # number of deprotonated groups
        if not no_deprotonated == None: 
            line+='%2d '%no_deprotonated
        else:
            line+=' - '
        # errors
        avr_error = 0.0
        for e in errors:
            line += '%8.5f '%e
            avr_error+=float(e)
        avr_error = avr_error/len(errors)
        
        line += '%8.5f '%avr_error
        line+='\n'
        # write output line
        self.report.write(line)
        self.report.flush()
        return
        





if __name__=='__main__':
    import sys
    if len(sys.argv)==2:
        W = wrapper(script=sys.argv[1])
    else:
        print 'Usage pKa_system_wrapper <script>'
        W = wrapper()



