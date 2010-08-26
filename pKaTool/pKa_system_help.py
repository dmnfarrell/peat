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

import sys

class system_help:

    def deactivate_callbacks(self):
        """Deactivate all callbacks from the scale widgets"""
        for group in self.groups:
            self.groups[group].intpka.configure(command=self.dummy)
            self.groups[group].active_btn.configure(command=self.dummy)
            self.groups[group].ab_button.configure(command=self.dummy)
            self.groups[group].style_button.configure(command=self.dummy)
            for group2 in self.groups[group].intenes.keys():
                self.groups[group].intenes[group2].configure(command=self.dummy)
        #
        # Unfolded form
        #
        if self.stability_var.get()=='on':
            for group in self.unfolded_groups:
                self.unfolded_groups[group].intpka.configure(command=self.dummy)
                self.unfolded_groups[group].active_btn.configure(command=self.dummy)
                self.unfolded_groups[group].ab_button.configure(command=self.dummy)
                self.unfolded_groups[group].style_button.configure(command=self.dummy)
                for group2 in self.unfolded_groups[group].intenes.keys():
                    self.unfolded_groups[group].intenes[group2].configure(command=self.dummy)
        return

    #
    # ----
    #
    
    def activate_callbacks(self):
        """Deactivate all callbacks from the scale widgets"""
        for group in self.groups:
            self.groups[group].intpka_scale.configure(command=self.update_scales)
            self.groups[group].active_btn.configure(command=self.update_pkasystem_curves)
            self.groups[group].ab_button.configure(command=self.update_pkasystem_curves)
            self.groups[group].style_button.configure(command=self.update_pkasystem_curves)
            for group2 in self.groups[group].intenes.keys():
                self.groups[group].intenes_scales[group2].configure(command=self.update_scales)
        #
        # Unfolded form
        #
        if self.stability_var.get()=='on':
            for group in self.unfolded_groups:
                self.unfolded_groups[group].intpka_scale.configure(command=self.update_scales)
                self.unfolded_groups[group].active_btn.configure(command=self.update_pkasystem_curves)
                self.unfolded_groups[group].ab_button.configure(command=self.update_pkasystem_curves)
                self.unfolded_groups[group].style_button.configure(command=self.update_pkasystem_curves)
                for group2 in self.unfolded_groups[group].intenes.keys():
                    self.unfolded_groups[group].intenes_scales[group2].configure(command=self.update_scales)
        return

#
# -----
#

class pKsensitivity:

    def change_dielectric(self):
        """Change the dielectric constant - i.e. scale all energies"""
        self.eps=Toplevel()
        self.eps.title('Change dielectric constant')
        self.eps.geometry('+200+200')
        Label(self.eps,text='Change dielectric constant by a factor').grid(row=0,column=0)
        self.eps_var=DoubleVar()
        self.eps_var.set(0.5)
        Entry(self.eps,textvariable=self.eps_var).grid(row=0,column=1)
        #
        # Check if we have model pKa values
        #
        self.has_modelpKs=1
        for group in self.groups.keys():
            if not self.groups[group].modelpK:
                self.has_modelpKs=None
                break
        if self.has_modelpKs:
            Label(self.eps,text='Changing intrinsic pKa values and interaction energies',bg='green').grid(row=1,column=0,columnspan=2)
        else:
            Label(self.eps,text='Cannot change intrinsic pKa values. Changing only interaction energies',bg='red').grid(row=1,column=0,columnspan=2)
        #
        # Submit button
        #
        go_button=Button(self.eps,text='Go',command=self.go_change_eps)
        go_button.grid(row=2,column=0,sticky='news',columnspan=2)
        return

    #
    # ----
    #
    def go_change_eps(self):
        """Scale all the pertubation energies"""
        self.stab_test_on=None
        factor=self.eps_var.get()
        #
        # Scale for all groups
        #
        for group in self.groups.keys():
            #
            # Explore the difference in the intrinsic pKa
            #
            if self.has_modelpKs:
                #
                # Vary as fraction of desolv+backgr (=intpka-modelpk)
                #
                this_var=(self.groups[group].intpka.get()-self.groups[group].modelpK)*factor
            else:
                this_var=0.0
                print 'Could not get model pKa value'
            #
            #
            old_intpka=self.groups[group].intpka.get()
            self.groups[group].intpka.set(old_intpka-float(this_var))
            #
            # Now change the interaction energies for this group
            #
            groups2=self.groups[group].intenes.keys()
            groups2.sort()
            for igroup in groups2:
                ene=self.groups[group].intenes[igroup].get()
                self.groups[group].intenes[igroup].set(factor*ene)

        #
        # normal mode
        #
        X=self.update_pkasystem_curves(doit=1)
        #
        # Should we update the display?
        #
        self.titwin.update()
        self.stab_test_on=1
        return
    #
    # ----
    #

    def sensitivity_test(self):
        """Vary system parameters and record the differences in the titration curves.
        This function opens the window for setting up the analysis"""
        #
        # Number of steps for varying intrinsic pKa values and interaction energies
        #
        intpka_steps=3
        intene_steps=3
        #
        # Open a window
        #
        self.Q=Toplevel()
        self.Q.title('System sensitivity analysis')
        import Window_geom
        Window_geom.set_geometry(self.titwin,self.Q)
        #self.Q.geometry('+200+200')
        #
        # Pulldown menu
        #
        self.menu=Menu(self.Q)
        #
        # File menu
        #
        self.monitor={}
        self.group_menu=Menu(self.menu,tearoff=0)
        groups=self.groups.keys()
        groups.sort()
        for group in groups:
            self.monitor[group]=IntVar()
            self.group_menu.add_checkbutton(label=group,
                                            variable=self.monitor[group],
                                            onvalue=1,offvalue=0)
        self.menu.add_cascade(label='Groups',menu=self.group_menu)
        self.Q.config(menu=self.menu)
        #
        # Write instructinos
        #
        row=0
        Label(self.Q,text='1. Use pulldown to select groups to monitor.').grid(row=row,column=0,sticky='news',columnspan=2)
        row=row+1
        Label(self.Q,text='2. Next select analysis method.').grid(row=row,column=0,sticky='news',columnspan=2)
        row=row+1
        #
        # Set the parameters
        #
        self.intene_var=DoubleVar()
        self.intene_var.set(20)
        self.intpka_var=DoubleVar()
        
        Label(self.Q,text='3. Specify variation parameters').grid(row=row,column=0,sticky='news',columnspan=2)
        row=row+1
        #
        # Test if we have model pKa values for all our groups. If we do, then we vary intrinsic pKa
        # values by a fraction by the desolv+backgr interaction energies. If not then we have no
        # choice but to vary the intrinsic pKa values by an absolute amount
        #
        self.has_modelpKs=1
        for group in groups:
            if not self.groups[group].modelpK:
                self.has_modelpKs=None
                break
        #
        # Well..
        #
        Label(self.Q,text='Vary intpKas by +/-').grid(row=row,column=0,sticky='news')
        if not self.has_modelpKs:
            self.intpka_var.set(0.5)
            Scale(self.Q,variable=self.intpka_var,resolution=0.1,from_=0.0,to=5.0,
                  label='pKa units',orient='horizontal').grid(row=row,column=1)
        else:
            self.intpka_var.set(50.0)
            Scale(self.Q,variable=self.intpka_var,resolution=5,from_=0.0,to=100.0,
                  label='%',orient='horizontal').grid(row=row,column=1)
        row=row+1
        #
        # Interaction energies are easier...
        #
        Label(self.Q,text='Vary interaction energies by +/-').grid(row=row,column=0,sticky='news')
        Scale(self.Q,variable=self.intene_var,resolution=5,from_=0.0,to=100,
              label='%',orient='horizontal').grid(row=row,column=1,sticky='news')
        row=row+1
        
        #
        # Full or limited
        #
        self.method=StringVar()
        self.method.set('single')
        Radiobutton(self.Q,variable=self.method,
                    text='Full combinatorial exploration',value='full').grid(row=row,column=0,sticky='news',columnspan=2)
        row=row+1
        Radiobutton(self.Q,variable=self.method,
                    text='Single energy at a time',value='single').grid(row=row,column=0,sticky='news',columnspan=2)
        row=row+1
        Label(self.Q,text='4. Press button to start analysis.').grid(row=row,column=0,sticky='news',columnspan=2)
        row=row+1
        #
        # Button for starting the analysis
        #
        self.go_button=Button(self.Q,text='Start analysis',command=self.go_sensitivity_test)
        self.go_button.grid(row=row,column=0,sticky='news',columnspan=2)
        #
        #
        #
        #
        return

    #
    # -----
    #

    def go_sensitivity_test(self):
        """Start the sensitivity test"""
        #
        # Did the user select any groups to monitor?
        #
        ok=None
        for id in self.monitor.keys():
            if self.monitor[id].get()==1:
                ok=1
                break
        if not ok:
            import tkMessageBox
            tkMessageBox.showwarning("No groups selected",
                                     'You have to select groups to monitor'
                                     ,parent=self.Q)
            return
        #
        # Delete the widgets from the setup window
        #
        self.menu.destroy()
        self.go_button.destroy()
        #
        # Put a wait label and a progress label
        #
        row=9
        wait_label=Label(self.Q,text='Sensitivity analysis running. Please wait',bg='green')
        wait_label.grid(row=row,column=0,columnspan=2,sticky='news')
        #
        row=row+1
        self.progress=StringVar()
        progress_label=Label(self.Q,textvariable=self.progress)
        progress_label.grid(row=row,column=0,sticky='news')
        #
        # Should we do the full analysis, or just the simple
        #
        self.draw=1
        if self.method.get()=='full':
            self.all_permutations()
        elif self.method.get()=='single':
            self.single_parameter_analysis()
        elif self.method.get()=='random':
            self.random_sensitivity_analysis()
        else:
            raise Exception('Incorrect value')
        #
        # OK, got all the results, now display them
        #
        wait_label.config(text='Displaying state and titration curve difference')
        yscrollbar=Scrollbar(self.Q,orient='vertical',width=10)
        yscrollbar.grid(row=row,column=1,sticky='nws')
        self.listbox=Listbox(self.Q,bg='white',
                             fg='black',
                             height=20,width=30,
                             yscrollcommand= yscrollbar.set,
                             selectmode=SINGLE)
        self.listbox.grid(row=row,column=0,sticky='news',columnspan=2)
        yscrollbar.config(command=self.listbox.yview)
        self.listbox.delete(0, END)
        count=1
        self.listbox.insert(END,'Reference state')
        for state,diff in self.curv_diffs[1:min(50,len(self.curv_diffs)-1)]:
            self.listbox.insert(END,'%d diff: %.2f' %(count,diff))
            count=count+1
        #
        # Back to the wt state
        #
        self.unpack_configuration(self.wt_state)
        #
        # bind the action to the display
        #
        self.listbox.bind('<ButtonRelease-1>',self.display_state)
        self.stab_test_on=None
        #
        # Quit Button
        #
        row=row+1
        Button(self.Q,text='Close window',command=self.close_Q).grid(row=row,column=0)
        return

    #
    # -----
    #

    def measure_change(self,reference=None,absolut=False):
        """Measure the difference in the titration curve as compared to the wt situration"""
        X=self.update_pkasystem_curves(doit=1)
        #
        # Should we update the display?
        #
        if self.draw:
            self.titwin.update()
        else:
            if self.counter%10==0:
                self.titwin.update()
        #
        # Get the reference
        #
        if not reference:
            reference=self.WT
        #
        this_diff=0.0
        self.counter=self.counter+1
        import time
        timenow=time.time()
        time_passed=timenow-self.starttime
        fraction_completed=float(self.counter)/float(self.do_states)
        time_to_finish=(1.0-fraction_completed)/(fraction_completed/time_passed)
        txt='%d seconds' %(time_to_finish)
        if time_to_finish/60.0/60.0/24.0>365:
            txt='%5.2f years' %(time_to_finish/60.0/60.0/24.0/365.0)  
        elif time_to_finish/60.0/60.0>24:
            txt='%5.2f days' %(time_to_finish/60.0/60.0/24.0)    
        elif time_to_finish/60.0>60:
            txt='%5.2f hours' %(time_to_finish/60.0/60.0)
        elif time_to_finish>60:
            txt='%5.2f minutes' %(time_to_finish/60.0)
        self.progress.set('%.2f %% done (%5.2e of %5.2e). Finish in %s' %(fraction_completed*100.0,self.counter,self.do_states,txt))
        #
        # Measure the change
        #
        for group_test in self.ids.keys():
            id=self.ids[group_test]
            if self.monitor[id].get()==1:
                for pH in X.pHvalues:
                    new_crg=X.prot_states[group_test][pH]
                    old_crg=reference.prot_states[group_test][pH]
                    if absolut:
                        this_diff=this_diff+abs(new_crg-old_crg)*self.pHstep.get()
                    else:
                        this_diff=this_diff+(new_crg-old_crg)*self.pHstep.get()
        #
        # Store this difference
        #
        if abs(this_diff)>0.0:
            state=self.pack_configuration()
            self.curv_diffs.append([state,this_diff])
            #
            # Sort by largest difference
            #
            
            self.curv_diffs.sort(lambda x, y: -cmp(abs(x[1]),abs(y[1])))
            if self.curv_diffs[-1][0]==self.wt_state:
                self.curv_diffs=[[self.wt_state,0.0]]+self.curv_diffs
                del self.curv_diffs[-1]
            else:
                self.curv_diffs=[[self.wt_state,0.0]]+self.curv_diffs
            #self.curv_diffs=self.curv_diffs[:100000]
        return this_diff

    #
    # -----
    #

    def single_parameter_analysis(self):
        """Vary every single parameter and record changes in the protonation states"""
        #
        # Store a copy of the wild type state
        #
        self.WT=self.update_pkasystem_curves()
        self.wt_state=self.pack_configuration()
        #
        # Start looping for all intrinsic pKa values
        #
        num_groups=len(self.groups.keys())
        self.tot_perms=2*num_groups+num_groups*(num_groups-1)
        self.do_states=self.tot_perms
        self.stab_test_on=1
        #
        self.curv_diffs=[]
        self.counter=0
        import time
        self.starttime=time.time()
        for group in self.groups.keys():
            #
            # Explore the difference in the intrinsic pKa
            #
            if self.has_modelpKs:
                #
                # Vary as fraction of desolv+backgr (=intpka-modelpk)
                #
                intpka_var=self.intpka_var.get()/100.0
                this_var=(self.groups[group].intpka.get()-self.groups[group].modelpK)*intpka_var
            else:
                this_var=self.intpka_var.get()
            #
            for diff in [-this_var,this_var]:
                self.unpack_configuration(self.wt_state) # Restore the wt state
                old_intpka=self.groups[group].intpka.get()
                self.groups[group].intpka.set(old_intpka+float(diff))
                self.measure_change()
            #
            # Now change the interaction energies for this group
            #
            intene_var=self.intene_var.get()/100.0
            groups2=self.groups[group].intenes.keys()
            groups2.sort()
            for igroup in groups2:
                for factor in [1.00-intene_var,1.00+intene_var]:
                    self.unpack_configuration(self.wt_state) # Restore the wt state
                    ene=self.groups[group].intenes[igroup].get()
                    self.groups[group].intenes[igroup].set(factor*ene)
                    self.measure_change()
        #
        # normal mode
        #
        self.stab_test_on=1
        return

    def random_sensitivity_analysis(self):
        """Randomly vary the setup"""
        #
        # Store a copy of the wild type state
        #
        self.WT=self.update_pkasystem_curves()
        self.wt_state=self.pack_configuration()
        #
        # Start looping for all intrinsic pKa values
        #
        num_groups=len(self.groups.keys())
        self.tot_perms=self.sens_steps
        self.do_states=self.tot_perms
        self.stab_test_on=1
        #
        self.curv_diffs=[]
        self.counter=0
        import time
        self.starttime=time.time()
        while self.counter<self.do_states:
            for group in self.groups.keys():
                #
                # Change intrinsic pKa
                #
                if self.has_modelpKs:
                    #
                    # Vary as fraction of desolv+backgr (=intpka-modelpk)
                    #
                    intpka_var=self.intpka_var.get()/100.0
                    this_var=(self.groups[group].intpka.get()-self.groups[group].modelpK)*intpka_var
                else:
                    this_var=self.intpka_var.get()
                #
                import random
                change=random.uniform(-this_var,this_var)
                #
                # Now change the interaction energies for this group
                #
                intene_var=self.intene_var.get()/100.0
                groups2=self.groups[group].intenes.keys()
                groups2.sort()
                for igroup in groups2:
                    change=random.uniform( 1.00-intene_var,1.00+intene_var)
                    ene=self.groups[group].intenes[igroup].get()
                    self.groups[group].intenes[igroup].set(change*ene)
            self.measure_change()
            self.unpack_configuration(self.wt_state) # Restore the wt state
        #
        # normal mode
        #
        self.stab_test_on=1
        return
  
                
    # ----
    #

    def all_permutations(self):
        #
        # Store a copy of all the intrinsic pKa values
        #
        intpka_backup={}
        for group in self.groups.keys():
            intpka_backup[group]=self.groups[group].intpka.get()
        #
        # Get the wild type titration curves
        #
        self.WT=self.update_pkasystem_curves()
        self.wt_state=self.pack_configuration()
        self.levels=0
        self.values=[]
        self.max_count=[]
        #
        # Vary all parameters systematically.
        # Intrinsic pka values +/- 0.5 
        # Interaction energies from +/- 20% 
        #
        self.numlevels=len(self.wt_state.split('_')[1:])
        self.stab_test_on=1
        for group in self.groups.keys():
            #
            # Intrinsic pKa value
            #
            values=[]
            if self.has_modelpKs:
                #
                # Vary as fraction of desolv+backgr (=intpka-modelpk)
                #
                intpka_var=self.intpka_var.get()/100.0
                this_var=(self.groups[group].intpka.get()-self.groups[group].modelpK)*intpka_var
            else:
                this_var=self.intpka_var.get()
            #
            for diff in [-this_var,this_var]:
                values.append(intpka_backup[group]+float(diff))
            self.max_count.append(len(values))
            self.values.append(values)
            self.levels=self.levels+1
            #
            # Interaction energies
            #
            intene_var=self.intene_var.get()/100.0
            groups2=self.groups[group].intenes.keys()
            groups2.sort()
            for igroup in groups2:
                values=[]
                for factor in [1.00-intene_var,1.00+intene_var]:
                    ene=self.groups[group].intenes[igroup].get()
                    values.append(ene*factor)
                #
                self.max_count.append(len(values))
                self.values.append(values)
                self.levels=self.levels+1
        #
        # How many permutations in total?
        #
        self.tot_perms=1
        for l in range(self.levels):
            self.tot_perms=self.tot_perms*self.max_count[l]
        #
        # Loop over all permutations
        #
        self.counter=0
        import time
        self.starttime=time.time()
        self.curv_diffs=[]
        #
        # Should we draw?
        #
        
        if self.tot_perms>100000 and not self.ignore_warnings:
            self.draw=None
            import tkMessageBox
            if not tkMessageBox.askyesno("Large system",
                                         "System is too large for sensitivity analysis\nContinue anyway?",
                                         parent=self.Q):
                self.Q.destroy()
                return
        #
        # ---
        #
        self.do_states=self.tot_perms
        if not hasattr(self,'sample'):
            self.sample=False
        if self.sample:
            #
            # Just probe random states to see what happens
            #
            self.do_states=self.sens_steps
            while self.counter<self.sens_steps:
                import random
                perm=random.randint(0,self.tot_perms-1)
                state=self.get_specs(perm)
                self.unpack_configuration(state)
                self.measure_change()
        else:
            while self.counter<min(self.tot_perms,100000):
                #
                # Get this configuration
                #
                perm=self.counter
                state=self.get_specs(perm)
                self.unpack_configuration(state)
                #
                # Calculate the titation curves and measure the change
                #
                self.measure_change()
        #
        self.stab_test_on=None
        return

    #
    # -----
    #

    def close_Q(self,event=None):
        """Close the Q window"""
        self.unpack_configuration(self.wt_state)
        self.update_pkasystem_curves()
        self.titwin.title('Titration curves [native state]')
        self.Q.destroy()
        return

    #
    # ------
    #

    def pack_configuration(self):
        """Pack the present configuration of the system into a string"""
        groups=self.groups.keys()
        groups.sort()
        s=''
        #
        # Intrinsic pKa values
        #
        for group in groups:
            s=s+'_%.2f' %(self.groups[group].intpka.get())
            #
            # Interaction energies
            #
            groups2=self.groups[group].intenes.keys()
            groups2.sort()
            for igroup in groups2:
                s=s+'_%.2f' %(self.groups[group].intenes[igroup].get())
        #
        # Done
        #
        return s

    #
    # ------
    #

    def unpack_configuration(self,s):
        """Given a state of the system in a string, set that state"""
        groups=self.groups.keys()
        groups.sort()
        s=s.split('_')
        s=s[1:]
        count=0
        #
        # Intrinsic pKa values
        #
        for group in groups:
            self.groups[group].intpka.set(float(s[count]))
            count=count+1
            #
            # Interaction energies
            #
            groups2=self.groups[group].intenes.keys()
            groups2.sort()
            for igroup in groups2:
                self.groups[group].intenes[igroup].set(float(s[count]))
                count=count+1
        return

    #
    # -----------
    #

    def get_specs(self,perm_number):
        """Given a permutation number get the configuration"""
        #
        base_number=[]
        base=1
        for level in range(self.levels-1,-1,-1):
            base_number.append(base)
            base=base*self.max_count[level]
        base_number.reverse()
        #
        # Find the configuration
        #
        conf=[]
        for base in base_number:
            frac=float(perm_number)/float(base)
            if frac>=1.0:
                num=int(frac)
                conf.append(num)
                perm_number=perm_number-num*base
            else:
                conf.append(0)
        #
        # Fill in the real values
        #
        configuration=''
        count=0
        for index in conf:
            value=self.values[count][index]
            import types
            if type(value) is types.IntType:
                value=float(value)/self.array_factor
            configuration=configuration+'_%.2f' %value
            count=count+1
        return configuration

    #
    # ----
    #

    def display_state(self,junk=None):
        """Display a particular state of the system"""
        for state in self.listbox.curselection():
            state_num=int(state)
            self.titwin.title('Titration curves [state: %d]' %state_num)
            state=self.curv_diffs[state_num][0]
            self.unpack_configuration(state)
        self.update_pkasystem_curves()
        self.update_scales_from_fit()
        return

        

#
# ----
#

class decompose(pKsensitivity):
    """Class for handling the decomposition of the system"""

  

    def decompose_system(self,event=None):
        """Switch off each group in turn, and record the change in a set of given titration curves"""
        #
        # Open a window
        #
        self.Q=Toplevel()
        self.Q.title('System decomposition')
        self.Q.geometry('+200+200')
        #
        # Just start since there are no parameters
        # 
        self.go_decompose()
        return

    #
    # -----
    #


    def go_decompose(self):
        """Actually do the decomposition analysis"""
        #
        # Put a wait label and a progress label
        #
        row=5
        wait_label=Label(self.Q,text='Decomposition analysis running. Please wait')
        wait_label.grid(row=row,column=0,sticky='news')
        #
        row=row+1
        self.progress=StringVar()
        progress_label=Label(self.Q,textvariable=self.progress)
        progress_label.grid(row=row,column=0,sticky='news')
        #
        #
        self.draw=1
        self.calc_decomposition()
        #
        #
        wait_label.destroy()
        progress_label.destroy()
        row=7
        #
        # Print it all
        #
        
        groups=self.ex_system.keys()
        groups.sort()
        
        
        #
        # print top labels
        #
        Label(self.Q,text='In_system energies').grid(row=row,column=1,columnspan=len(groups))
        row=row+1
        row=self.print_matrix(row,self.in_system,groups)+1
        Label(self.Q,text='-----------------------------------------------').grid(row=row,column=0,columnspan=len(groups))
        row=row+1
        #
        #
        Label(self.Q,text='Ex_system energies').grid(row=row,column=1,columnspan=len(groups))
        row=row+1
        row=self.print_matrix(row,self.ex_system,groups)+1
        Label(self.Q,text='-----------------------------------------------').grid(row=row,column=0,columnspan=len(groups))
        row=row+1
        #
        # Calculate the difference (in_system-ex_system)
        #
        diff={}
        for group1 in groups:
            if not diff.has_key(group1):
                diff[group1]={}
            for group2 in groups:
                if group2!=group1:
                    diffs=self.in_system[group1][group2]-self.ex_system[group1][group2]
                    if abs(self.ex_system[group1][group2])!=0.0:
                        sign=self.ex_system[group1][group2]/abs(self.ex_system[group1][group2])
                        diff[group1][group2]=diffs/sign
                    else:
                        diff[group1][group2]='N/A'
        #
        Label(self.Q,text='Diffs').grid(row=row,column=1,columnspan=len(groups))
        row=row+1
        self.print_matrix(row,diff,groups)
        return
    #
    # ----
    #

    def print_matrix(self,row,array,groups):
        column_num=2
        g_count=0
        Label(self.Q,text='Removed group').grid(row=row,column=1,columnspan=len(groups))
        row=row+1
        for monitor_group in groups:
            colour=self.colour_order[g_count%len(self.colour_order)]
            if colour=='black':
                fgc='white'
            else:
                fgc='black'
            Label(self.Q,text=monitor_group,bg=colour,fg=fgc).grid(row=row,column=column_num,sticky='news')
            column_num=column_num+1
            g_count=g_count+1
        #
        # ----
        #
        row=row+1
        #
        # Label
        #
        Label(self.Q,text='Monitored group').grid(row=row,column=0,rowspan=len(groups))
        g_count=0
        for monitor_group in groups:
            column=1
            #
            # Get the colour of the label right
            #
            colour=self.colour_order[g_count%len(self.colour_order)]
            if colour=='black':
                fgc='white'
            else:
                fgc='black'
            Label(self.Q,text=monitor_group,bg=colour,fg=fgc).grid(row=row,column=column,sticky='news')
            g_count=g_count+1
            column=column+1
            #
            # Show the results
            #
            for group2 in groups:
                if array[monitor_group].has_key(group2):
                    if array[monitor_group][group2]!='N/A':
                        Label(self.Q,text='%.2f' %array[monitor_group][group2]).grid(row=row,column=column,sticky='news')
                    else:
                        Label(self.Q,text='%s' %array[monitor_group][group2]).grid(row=row,column=column,sticky='news')
                else:
                    Label(self.Q,text='-').grid(row=row,column=column,sticky='news')
                column=column+1
            row=row+1
        return row

    #
    # ----
    #


    def calc_decomposition(self):
        """Run the insystem and exsystem decomposition energies"""

        #
        # Store a copy of the wild type state
        #
        self.WT=self.update_pkasystem_curves()
        self.wt_state=self.pack_configuration()
        #
        # Start looping over all groups
        #
        num_groups=len(self.groups.keys())
        self.tot_perms=2*num_groups+num_groups*(num_groups-1)
        #
        self.curv_diffs=[]
        self.counter=0
        self.in_system={}
        self.ex_system={}
        groups=self.groups.keys()
        groups.sort()
        for monitor_group in groups:
            #
            # Make sure we get the correct titration curve difference
            #
            self.monitor={}
            for group in groups:
                self.monitor[group]=IntVar()
                self.monitor[group].set(0)
                if group==monitor_group:
                    self.monitor[group].set(1)
            #
            # Swith off all other groups
            #
            # Restore the wt state
            for group_t in groups:
                self.groups[group_t].active.set(1)
            #
            for group2 in groups:
                if monitor_group!=group2:
                    self.groups[group2].active.set(0)
            alone_tit=self.update_pkasystem_curves(doit=1)
            #
            # Prepare the dictionaries
            #
            if not self.in_system.has_key(monitor_group):
                self.in_system[monitor_group]={}
            if not self.ex_system.has_key(monitor_group):
                self.ex_system[monitor_group]={}
            #
            # Loop over the second groups for getting the ex_system titrations
            #
            for group_off in groups:
                if group_off!=monitor_group:
                    #
                    # Turn on this group
                    #
                    self.groups[group_off].active.set(1)
                    #
                    # Get the difference between the isolated titration curve
                    # and the system consisting of only two groups
                    #
                    ex_sys_diff=-self.measure_change(alone_tit,absolut=None)
                    self.ex_system[monitor_group][group_off]=ex_sys_diff
                    #
                    # Turn the group off again
                    #
                    self.groups[group_off].active.set(0)
            #
            # Go on to calculating the in_system energies
            #
            # Restore the wt state
            for group_t in groups:
                self.groups[group_t].active.set(1)
            #
            for group_off in groups:
                if group_off!=monitor_group:
                    #
                    # Turn off one group
                    #
                    self.groups[group_off].active.set(0)
                    diff=self.measure_change(absolut=None)
                    self.in_system[monitor_group][group_off]=diff
                    #
                    # Turn the group back on
                    #
                    self.groups[group_off].active.set(1)
        #
        # Get the wild type state back
        #
        # Restore the wt state
        for group_t in groups:
            self.groups[group_t].active.set(1)
        return
    
