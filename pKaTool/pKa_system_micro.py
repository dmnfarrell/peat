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
from numpy import *
from numpy.linalg import * 
inverse=inv
Float=float
    
import math

from titration_class import *
import pKa_system_micro2


class Micro_states(pKa_system_micro2.Micro_states_helper):

    def update_microstates(self,X):
        """Update the curves describing the microscopic states"""
        self.X = X
        phvals=X.all_states.keys()
        self.init_kcat_array()

        if getattr(self,'activity_data',None):
            for ph in self.activity_data:
                if not ph in phvals:
                    phvals.append(ph)

        phvals.sort()
        self.phvals=phvals
        #
        # Delete all lines
        #
        for obj in self.state_lines.keys():
            self.state_win_cv.delete(obj)
            del self.state_lines[obj]
        #
        # Get the states
        #
        self.states=states=X.all_states[phvals[0]].keys()
        states.sort()
        #
        # Do we have the window open?
        #
        if not getattr(self,'state_win',None):
            self.init_win(X)
        #
        # Draw the states
        #
        state_count=0

        for state in states:
            #
            # Get the colour
            #
            colour=self.colour_order[state_count%len(self.colour_order)]
            fg_col='black'
            if colour=='black' or colour=='blue':
                fg_col='white'
            if colour=='green':
                colour='cyan'
            state_count=state_count+1
            #
            # Draw the state
            #
            lastpH=phvals[0]
            lastcrg=X.all_states[lastpH][state]['pop']
            #
            # Normalise the population
            #
            pops=[]
            for pH in phvals[1:]:
                pops.append(X.all_states[pH][state]['pop'])
            max_pop=max(pops)
            if max_pop>0.0:
                factor=1.0/max_pop
            else:
                factor=1.0
            self.state_pop[state].set('Max pop: %5.2e' %max_pop)
            #
            # Draw the state if the user wants to
            #
            charges=[]
            for pH in phvals[1:]:
                lastx,lasty=self.get_xy(lastpH,lastcrg)
                crg=X.all_states[pH][state]['pop']
                x,y=self.get_xy(pH,crg)
                if self.state_box[state].get()==1:
                    self.state_lines[(self.state_win_cv.create_line(lastx,lasty,float(x),float(y),
                                                                    fill=colour,
                                                                    width=self.linewidth))]=1
                #
                # Scale the maximum population to 1.0 for pKa fitting
                #
                charges.append([pH,crg*factor])
                lastcrg=crg
                lastpH=pH
            #
            # Determine the pK1/2 values for this state
            #
            if self.determine_pkhalf_checkvar.get()==1:
                pkas=self.determine_pkhalf(charges,1)
                if len(pkas)>0:
                    self.pk1[state].set('%5.2f' %pkas[0])
                    if len(pkas)>1:
                        self.pk2[state].set('%5.2f' %pkas[1])
            else:
                self.pk1[state].set('NA')
                self.pk2[state].set('NA')
        #
        # Draw the primary CCPS
        #
        act_prof = self.get_act_prof_and_tc_and_ftir()[0]
        if act_prof!={}:
            lastpH=phvals[0]
            last_activity=act_prof[lastpH]
            for pH in phvals[1:]:
                lastx,lasty=self.get_xy(lastpH,last_activity)
                activity=act_prof[pH]
                x,y=self.get_xy(pH,activity)
                self.state_lines[(self.state_win_cv.create_line(lastx,lasty,float(x),float(y),
                                                                fill='green',
                                                                width=3))]=1
                last_activity=activity
                lastpH=pH
        #
        # Secondary CCPS
        #
        act_prof = self.get_act_prof_and_tc_and_ftir()[-1]
        if act_prof!={}:
            lastpH=phvals[0]
            last_activity=act_prof[lastpH]
            for pH in phvals[1:]:
                lastx,lasty=self.get_xy(lastpH,last_activity)
                activity=act_prof[pH]
                x,y=self.get_xy(pH,activity)
                self.state_lines[(self.state_win_cv.create_line(lastx,lasty,float(x),float(y),
                                                                fill='orange',
                                                                width=3))]=1
                last_activity=activity
                lastpH=pH
        #
        # Draw the primary experimental pH-activity profile
        #
        if getattr(self,'activity_data',None):    
            phs = self.activity_data.keys()
            phs.sort()
            for ph in phs:
                activity = self.activity_data[ph]
                x,y=self.get_xy(float(ph),float(activity))
                handle=self.state_win_cv.create_oval(x-2,y-2,x+2,y+2,fill='green')
                self.state_lines[handle]=1
        #
        # Draw the secondary experimental pH-activity profile
        #
        if getattr(self,'secondary_activity_data',None):    
            phs = self.secondary_activity_data.keys()
            phs.sort()
            for ph in phs:
                activity = self.secondary_activity_data[ph]
                x,y=self.get_xy(float(ph),float(activity))
                handle=self.state_win_cv.create_oval(x-2,y-2,x+2,y+2,fill='orange')
                self.state_lines[handle]=1

                    
        return

    #
    # ----
    #
     
    def determine_pkhalf(self,charges,acid_base):
        """from a list of [pH,charge] determine two pKa values that define the population"""
        import copy
        #
        # First find out if there are one or two pKa values
        #
        pkas=self.determine_pkhalf_simple(charges,acid_base)
        if len(pkas)==2:
            #
            # Function for two pKa values
            #
            self.funktion="(1.0/(math.pow(10,-pH)/math.pow(10,-pKa1)+1+math.pow(10,-pKa2)/math.pow(10,-pH)))"
            self.funktion=self.funktion.replace('pKa1','self.pkvars[0][1]')
            self.funktion=self.funktion.replace('pKa2','self.pkvars[1][1]')
            self.pkvars=[[4.0,4.0],[6.0,6.0]]
        else:
            if charges[0][1]>0.2:
                
                self.funktion="(1.0-(1.0/(1.0+math.pow(10,pKa1-pH))))"
            else:
                self.funktion="(1.0/(1.0+math.pow(10,pKa1-pH)))"
            self.funktion=self.funktion.replace('pKa1','self.pkvars[0][1]')
            self.pkvars=[[4.0,4.0]]
        #
        # Fit the function
        #
        old_diff=10000
        for x in range(200):
            self.fit_LM_micro(charges)
            #
            # Check convergence
            # 
            now_diff=self.calc_sq_diff_micro(charges)
            if (abs(now_diff-old_diff)/old_diff)<0.0001:
                break
            else:
                old_diff=now_diff
        #
        # Compile the results
        #
        rets=[]
        for pka in self.pkvars:
            rets.append(pka[1])

        return rets

    #
    # ----
    #

    def determine_pkhalf_simple(self,charges,acid_base):
        """from a list of [pH,charge] determine the pK1/2"""
        pka=[]
        last_crg=charges[0][1]
        phstep=float(charges[1][0]-charges[0][0])
        for ph,crg in charges:
            if acid_base==1:
                if crg<=0.5 and last_crg>0.5:
                    pka.append((last_crg-0.5)/(last_crg-crg)*phstep+(ph-phstep))
                elif crg>=0.5 and last_crg<0.5:
                    pka.append((last_crg-0.5)/(last_crg-crg)*phstep+(ph-phstep))
            else:
                if crg<=-0.5 and last_crg>-0.5:
                    pka.append((last_crg-(-0.5))/(last_crg-crg)*phstep+(ph-phstep))
            last_crg=crg
        return pka


    #
    # ----
    #

    def fit_to_ph_activity_profile(self,silent=0):

        if not getattr(self,'states',None):
            self.micro_var.set(1.0)
            self.update_pkasystem_curves()

        if not getattr(self,'activity_data',None):
            import tkMessageBox
            tkMessageBox.showwarning('No pH-activity profile loaded',
                                     'Load pH-activity profile first')
            return


        no_active = 0
        for state in self.states:
            if self.act_state[state].get()==1:
                no_active = no_active+1

        if no_active == 0:
            import tkMessageBox
            tkMessageBox.showwarning('No active states selected',
                                     'Select active micro states first')
            return

        if self.include_tcs_in_fit == 1:
            if not getattr(self,'titration_data',None):
                import tkMessageBox
                tkMessageBox.showwarning('No titration data loaded',
                                     'Load titration data first')
                return
            else:
                self.exp_data=titration_curve(self.titration_data)
        
        if self.include_ftir_in_fit == 1:
            if not getattr(self,'FTIR_win',None):
                import tkMessageBox
                tkMessageBox.showwarning('No FTIR data loaded',
                                     'Load FTIR data first')
                return    
        #
        # Make copies of all the variables
        #
        self.vars=[]
        for group in self.groups:
            self.vars.append(self.groups[group].intpka)
            for group2 in self.groups[group].intenes.keys():
                self.vars.append(self.groups[group].intenes[group2])
        
        if self.include_ftir_in_fit == 1:
            self.vars.append(self.FTIR_win.c)
            self.vars.append(self.FTIR_win.offset)
        
        
        # Set damper
        #
        self.LM_damper = DoubleVar()
        self.LM_damper.set(0.1)
        #
        # Open the window
        #
        self.keep_running_activity=1
        if not silent:
            self.count_win=Toplevel()
            self.count_win.geometry('100x100+200+200')
            self.count_win.title('Fit in progress')
            self.count=IntVar()
            self.count.set(0)
            Label(self.count_win,text='Iteration #').grid(row=0,column=0) 
            Label(self.count_win,textvariable=self.count).grid(row=0,column=1)
            Button(self.count_win,text='Stop fit',command=self.stop_fit_activity).grid(row=2,column=0,columnspan=2)
            self.count_win.update_idletasks()



        old_diff=0.0000001
        print 'Step    Diff   LM_damper'
        now_diff=self.activity_diff()
        print '%5d  %6.4f  %6.4f' %(0, now_diff, self.LM_damper.get())
        for x in range(1,100):###100
            if self.keep_running_activity == 0:
                break
            self.fit_LM_activity()
            #
            # Update more stuff
            #
            if not silent:
                self.count.set(self.count.get()+1)
                self.count_win.update()
                self.update_scales_from_fit()
                self.titwin.update()
                if getattr(self,'FTIR_win',None):
                    self.FTIR_win.draw_fit()

            now_diff=self.activity_diff()
            #
            # Check convergence
            #
            print '%5d  %6.4f  %6.4f' %(x, now_diff, self.LM_damper.get())

            if abs(now_diff-old_diff)<0.00005:
                print 'Converged',now_diff
                break
            else:
                old_diff=now_diff

            if not silent:
                self.count.set(self.count.get()+1)
                self.count_win.update()
            self.update_scales_from_fit()
            self.titwin.update()
            if getattr(self,'FTIR_win',None):
                self.FTIR_win.draw_fit()

            now_diff=self.activity_diff()
           

        if not silent:
            self.count_win.destroy()
            self.update_pkasystem_curves()

        return



    def fit_LM_activity(self):
        """Do Levenberg-Marquardt fitting"""
        J,E =self.get_jacobian_activity()
        JT = transpose(J)
        JTE = dot(JT,E)
        JTJ = dot(JT,J)
        JTJd = JTJ + self.LM_damper.get()*identity(shape(JTJ)[0])
        invJTJd = inv(JTJd)
        q = -dot(JTE,invJTJd)

#        print 'q',q
#        print 'J',J
#        print 'E',E
        out1 ='bv '
        out2 ='av '
        for var in range(len(self.vars)):
            out1 += '%4.6f '%self.vars[var].get()
            self.vars[var].set(self.vars[var].get()+q[var])
            out2 += '%4.6f '%self.vars[var].get()
#        print out1
#        print out2
        return



    def get_jacobian_activity(self):
        """Get the Jacobian matrix and errors of the data points"""
        #
        # Get the number of data points
        #
        no_data_points = len(self.activity_data.keys())
        if self.include_tcs_in_fit == 1:
            for group in self.exp_data.curves.keys():
                no_data_points = no_data_points + len(self.exp_data.curves[group].keys())
       
        if self.include_ftir_in_fit == 1:
            no_data_points = no_data_points + len(self.FTIR_win.ftir_data.keys())
        #
        #
        #
        errors = resize(array(0,float),[no_data_points])        
        jacobian = resize(array(0,float),[no_data_points,len(self.vars)])
        #
        # Precalculate the variation of all parameters
        #
        now=self.get_act_prof_and_tc_and_ftir()
        variations=[]
        step = 1e-8

        for var in range(len(self.vars)):
            self.vars[var].set(self.vars[var].get()+step)
            variations.append(self.get_act_prof_and_tc_and_ftir())
            self.vars[var].set(self.vars[var].get()-step)
        #
        # construct jacobian
        #
        data_id=0

        #
        # activity data
        #
        for ph in self.activity_data.keys():
            data_point = self.activity_data[ph]
            x=float(ph)
            y=float(data_point)
                
            errors[data_id] = y-now[0][ph]
            #
            # Find the derivative at this ph (x) for this data point
            #
            diff=resize(array(0,float),[len(self.vars)])
            count=0
            for variation in variations:
                diff[count]=(now[0][ph]-variation[0][ph])/step
                count=count+1
            jacobian[data_id]=diff
            data_id=data_id+1
        #
        # titration data
        #
        if self.include_tcs_in_fit == 1:
            for group in self.exp_data.curves.keys():
                for ph in self.exp_data.curves[group].keys():
                    data_point = self.exp_data.curves[group][ph]
                    x=float(ph)
                    y=float(data_point)
                
                    errors[data_id] = y-self.titration_curves[group][ph]
                    #
                    # Find the derivative at this ph (x) for this data point
                    #
                    diff=resize(array(0,float),[len(self.vars)])
                    count=0
                    for variation in variations:
                        diff[count]=(now[1].curves[group][ph]-variation[1].curves[group][ph])/step
                        count=count+1
                    jacobian[data_id]=diff
                    data_id=data_id+1
        #
        # FTIR data
        #
        if self.include_ftir_in_fit == 1:
            sum_of_tits = self.FTIR_win.get_sum_of_tits()
            for ph in self.FTIR_win.ftir_data.keys():
                errors[data_id] = self.FTIR_win.ftir_data[ph]-sum_of_tits[ph]            
                #
                # Find the derivative at this ph (x) for this data point
                #
                diff=resize(array(0,float),[len(self.vars)])
                count=0
                for variation in variations:
                    #print '%f - %f = %f:' %(now[ph],variation[ph],now[ph]-variation[ph])
                    diff[count]=(now[2][ph]-variation[2][ph])/step
                    count=count+1
                jacobian[data_id]=diff
                data_id=data_id+1
                        
                    
        return jacobian,errors


    #
    # ---
    #
    def fit_to_curves_and_ph_activity_profile(self):
        self.include_tcs_in_fit = 1
        self.fit_to_ph_activity_profile()
        self.include_tcs_in_fit = 0
        return

        
    def fit_to_ftir_and_ph_activity_profile(self):
        self.include_ftir_in_fit = 1
        self.fit_to_ph_activity_profile()
        self.include_ftir_in_fit = 0
        return

    #
    # ----
    #

    def calc_sq_diff_micro(self,charges):
        """Calculate the square difference"""
        diff=0.0
        import math
        for pH,crg in charges:
            y=crg
            fit_val=eval(self.funktion)
            diff=diff+math.pow(fit_val-y,2)
        return diff

    #
    # ------
    #

    def fit_LM_micro(self,charges):
        """Do Levenberg-Marquardt fitting"""
        self.LM_damper=1.0
        J,E =self.get_jacobian_micro(charges)
        JT = transpose(J)
        JTE = dot(JT,E)
        JTJ = dot(JT,J)
        JTJd = JTJ + self.LM_damper*identity(shape(JTJ)[0])
        invJTJd = inverse(JTJd)
        q = -dot(JTE,invJTJd)
        #
        # Update the variables
        #
        for var in range(len(self.pkvars)):
            self.pkvars[var][1]=self.pkvars[var][1]+q[var]
        return

    #
    # -----
    #

    def get_jacobian_micro(self,charges):
        """Get the Jacobian matrix and errors of the data points"""
        no_data_points = len(charges)
        
        errors = resize(array(0,float),[no_data_points])        
        jacobian = resize(array(0,float),[no_data_points,len(self.pkvars)])
        temp = resize(array(0,float),[len(self.pkvars)])
        #
        #
        #
        count=0
        for pH,crg in charges:
            errors[count] = crg-eval(self.funktion)
            #
            #
            #
            for var in range(len(self.pkvars)):
                now=eval(self.funktion)
                self.pkvars[var][1]=self.pkvars[var][1]+0.1
                temp[var]=(now-eval(self.funktion))/0.1
                self.pkvars[var][1]=self.pkvars[var][1]-0.1
            #
            #
            #
            jacobian[count]=temp
            #
            # Update the counter
            #
            count=count+1
        return jacobian,errors

    #
    # ----
    #
    
    def init_states(self):
        """Get the states"""
        self.phvals=self.pKa_calc_instance.all_states.keys()
        self.phvals.sort()
        first_pH=self.phvals[0]
        self.states=self.pKa_calc_instance.all_states[first_pH].keys()
        self.states.sort()
        self.last_pH=self.phvals[-2]
        return
    
    #
    # ----
    #

    def init_win(self,pKa_calc_instance):
        """Initialise the windows for displaying and manipulating microscopic states"""
        #
        # Get a desription of all states
        #
        self.init_states()
        self.X=pKa_calc_instance
        #self.lastpH=lastpH
        self.include_tcs_in_fit = 0
        self.include_ftir_in_fit = 0
        #
        # Open the window
        #
        self.state_win=Toplevel()
        self.state_win.title('Microscopic states')
        self.state_win.geometry('+20+300')

        #
        # Checkboxes
        #
        self.act_state={}
        self.act_state2={}
        self.state_box={}
        self.state_pop={}
        self.pk1={}
        self.pk2={}
        state_count=0
        #
        # Place labels, and action buttons
        #
        Label(self.state_win,text='Show').grid(row=2,column=0)
        Label(self.state_win,text='State description (1: charged)').grid(row=2,column=1,columnspan=2)
        Label(self.state_win,text='pK1').grid(row=2,column=3,sticky='news')
        Label(self.state_win,text='pK2').grid(row=2,column=4,sticky='news')
        Label(self.state_win,text='CCPS').grid(row=2,column=5,sticky='news')
        Label(self.state_win,text="2nd").grid(row=2,column=6,sticky='news')
        Label(self.state_win,text='pH-dependent population of states').grid(row=1,column=7,sticky='news')
        Label(self.state_win,text='pH-activity profile as defined by Active boxes is shown in green').grid(row=2,column=7,sticky='news')
        offset=3
        #
        # Open the canvas
        #
        self.state_win_cv=Canvas(self.state_win,
                                 width=900,
                                 height=250,
                                 bg='white')
        self.state_win_cv.grid(row=offset,column=7,rowspan=10)
        #
        # Draw the axes
        #
        # pH axis
        self.micro_startx,self.micro_starty=self.get_xy(0.0,0.0)
        self.micro_endx=860
        self.micro_endy=40
        self.pH_axis(self.state_win_cv,
                     self.micro_startx,
                     self.micro_starty,
                     self.micro_endx,
                     self.micro_endy,
                     'pH')
        #
        # population axis
        #
        self.state_win_cv.create_line(self.micro_startx,
                                      self.micro_starty,
                                      self.micro_startx,
                                      self.micro_endy,fill='black',width=self.linewidth)
        self.state_win_cv.create_text(self.micro_startx,self.micro_endy-15,text='Population',fill='black')
        #
        # Tick marks and tick labels
        #
        for tickcrg in range(0,110,10):
            dummy,y=self.get_xy(0.0,float(tickcrg)/100.0)
            self.state_win_cv.create_line(self.micro_startx,y,self.micro_startx-5,y,fill='black',width=self.linewidth)
            self.state_win_cv.create_text(self.micro_startx-25,y,text='%5.2f' %(float(tickcrg)/100.0),fill='black')
        #
        # Write each state and the corresponding checkboxes
        #
        for state in self.states:
            #
            # Get colour
            #
            colour=self.colour_order[state_count%len(self.colour_order)]
            fg_col='black'
            if colour=='black' or colour=='blue':
                fg_col='white'
            if colour=='green':
                colour='cyan'
            state_count=state_count+1
            #
            # The box
            #
            self.state_box[state]=IntVar()
            self.state_box[state].set(1)
            Checkbutton(self.state_win,variable=self.state_box[state],
                        onvalue=1,offvalue=0,command=self.update_pkasystem_curves).grid(row=state+offset,column=0)
            #
            # Write the definition of the state
            #
            obj=Label(self.state_win,
                      text="state %d, def: %s" %(state,self.X.all_states[self.last_pH][state]['def']),
                      bg=colour,fg=fg_col)
            obj.grid(row=state+offset,column=1)
            #
            # Write the population
            #
            self.state_pop[state]=StringVar()
            self.state_pop[state].set('Uninitialised')
            obj=Label(self.state_win,
                      textvariable=self.state_pop[state],
                      bg=colour,fg=fg_col)
            obj.grid(row=state+offset,column=2)
            #
            # pK1 and pK2 values
            #
            self.pk1[state]=StringVar()
            self.pk1[state].set('NA')
            obj=Label(self.state_win,
                      textvariable=self.pk1[state])
            obj.grid(row=state+offset,column=3)
            # pK2
            self.pk2[state]=StringVar()
            self.pk2[state].set('NA')
            obj=Label(self.state_win,
                      textvariable=self.pk2[state])
            obj.grid(row=state+offset,column=4)
            #
            # The checkboxes for defining the active states
            #
            self.act_state[state]=IntVar()
            self.act_state[state].set(0)
            Checkbutton(self.state_win,variable=self.act_state[state],
                        onvalue=1,offvalue=0,command=self.update_pkasystem_curves).grid(row=state+offset,column=5)
            #
            # Secondary CCPS
            #
            self.act_state2[state]=IntVar()
            self.act_state2[state].set(0)
            Checkbutton(self.state_win,variable=self.act_state2[state],
                        onvalue=1,offvalue=0,command=self.update_pkasystem_curves).grid(row=state+offset,column=6)
        #
        # Write a checkbox for calculating the pKa values for the microstates
        #
        Label(self.state_win,text='Find pKa values for microstates:').grid(row=len(self.states)+offset,column=0,columnspan=3)
        self.determine_pkhalf_checkvar=IntVar()
        self.determine_pkhalf_checkvar.set(0)
        Checkbutton(self.state_win,variable=self.determine_pkhalf_checkvar,
                    onvalue=1,offvalue=0,command=self.update_pkasystem_curves).grid(row=len(self.states)+offset,column=3)
        Label(self.state_win,text='Disable pKa value fitting to speed up calculations').grid(row=len(self.states)+offset+1,
                                                                                             column=0,columnspan=4)
        #
        # Button for transferring to Ekin
        #
        Button(self.state_win,text='Copy protonation state population to Ekin',command=self.copy_to_Ekin_micro).grid(row=len(self.states)+offset+2,column=0,columnspan=4,sticky='news')
        Button(self.state_win,text='Evaluate fit of pH activity profile',command=self.evaluate_fit_ph_activity).grid(row=len(self.states)+offset+3,column=0,columnspan=4,sticky='news')
        Button(self.state_win,text='Dump graphs to file',command=self.print2file_micro).grid(row=len(self.states)+offset+4,column=0,columnspan=4,sticky='news')
        Button(self.state_win,text='Save CCPS curve to file',command=self.saveCCPS).grid(row=len(self.states)+offset+5,column=0,columnspan=4,sticky='news')
        #
        # Make sure we close nicely
        #
        self.state_win.protocol("WM_DELETE_WINDOW",self.close_state_win)
        return

    #
    # -----
    #

    def close_state_win(self):
        """Set the vars so we now the state_win is closed"""
        if getattr(self,'state_win',None):
            self.state_win.destroy()
            self.state_win=None
            self.micro_var.set(0)
            self.state_lines={}
        return

    #
    # -----
    #

    def saveCCPS(self):
        """Save the CCPS curve to file"""
        act_prof,curve,ftir,act_prof2=self.get_act_prof_and_tc_and_ftir()
        import sys, os, tkFileDialog
        if not self.lastdir:
            self.lastdir=os.getcwd()
        filename=tkFileDialog.asksaveasfilename(defaultextension='.txt',
                                                initialdir=self.lastdir,
                                                filetypes=[("Text file","*.txt"),("All files","*.*")],
                                                parent=self.state_win)
        if filename:
            fd=open(filename,'w')
            pHs=act_prof.keys()
            pHs.sort()
            fd.write('pH\tActivity\n')
            for pH in pHs:
                fd.write('%7.3f\t%7.3f\n' %(pH,act_prof[pH]))
            fd.close()  
        return

    #
    # -----
    #

    def print2file_micro(self):
        #
        # Print Canvas to file
        #
        import sys, os, tkFileDialog
        if not self.lastdir:
            self.lastdir=os.getcwd()
        filename=tkFileDialog.asksaveasfilename(defaultextension='.ps',
                                                initialdir=self.lastdir,
                                                filetypes=[("Postscript files","*.ps"),("All files","*.*")],
                                                parent=self.state_win)
        if filename:
            self.write_psfile_micro(filename)
        else:
            return
        return
    
    #
    # --------------------
    #    
   
    def write_psfile_micro(self,filename):
        #
        # Dump the Canvas
        #
        import os
        self.lastdir=os.path.split(filename)[0]
        if filename[-3:]!='.ps':
            filename=filename+'.ps'
        self.state_win_cv.postscript(colormode='color',file=filename)
        return

    #
    # ----
    #

    def copy_to_Ekin_micro(self):
        """Copy a population curve to the Ekin facility of EAT_DB"""
        import os,sys
        import EAT_DB.Ekin
        #
        # Pick a group
        #
        self.pick_group_micro=Toplevel()
        self.pick_group_micro.title('Pick a protonation state')
        self.pick_group_micro.geometry('+200+200')
        self.state_picked=IntVar()
        count=0
        for state in self.states:
            #
            # Get colour
            #
            colour=self.colour_order[count%len(self.colour_order)]
            fg_col='black'
            if colour=='black' or colour=='blue':
                fg_col='white'
            if colour=='green':
                colour='cyan'
            #
            # Put the radiobutton in the window
            #
            Radiobutton(self.pick_group_micro,text="state %d, def: %s" %(state,self.X.all_states[self.lastpH][state]['def']),
                        variable=self.state_picked,
                        value=count).grid(row=count,column=0)
            count=count+1
        self.state_picked.set(0)
        Button(self.pick_group_micro,text='Copy population',command=self.copy_group_micro).grid(row=count,column=0)
        Button(self.pick_group_micro,text='Cancel',command=self.cancel_copy_group_micro).grid(row=count,column=1)
        return

    #
    # ----
    #
    
    def copy_group_micro(self,event=None):
        """Get the population curve and send it to Ekin"""
        #
        # Get the data and reformat it
        #
        state=self.state_picked.get()
        data={}
        for pH in self.X_states.all_states.keys():
            data[pH]=self.X_states.all_states[pH][state]['pop']
        new_data={}
        new_data[0]={}
        new_data[1]={}
        count=0
        pHs=data.keys()
        pHs.sort()
        for pH in pHs:
            new_data[0][count]=pH
            new_data[1][count]=data[pH]
            count=count+1
        #
        # Open Ekin, and load the data
        #
        import os,sys
        import EAT_DB.Ekin
        EK=EAT_DB.Ekin.Ekin(parent=self)
        EK.pass_data(new_data,'pKaSystem population')
        #
        # Destroy the little window
        #
        self.pick_group_micro.destroy()
        return

    #
    # ----
    #

    def get_CCPSs(self):
        """Return a list of the protonation states chosen as being the CCPS"""
        X,pKa_values,prot_states=self.calc_pKas_from_scales(self.groups)
        curve=titration_curve(X.prot_states)
        chosen_states=[]
        for state in self.states:
            if self.act_state[state].get()==1:
                chosen_states.append(state)
        return chosen_states

    #
    # ----
    #

    def get_act_prof_and_tc_and_ftir(self):
        """Get the pH-activity profile and titration curves"""
        #
        # titration data
        #
        X,pKa_values,prot_states=self.calc_pKas_from_scales(self.groups)
        curve=titration_curve(X.prot_states)
        #
        # Primary activity data
        #
        act_prof1={}
        states=0
        for pH in self.phvals:
            act=0.0
            for state in self.states:
                if self.act_state[state].get()==1:
                    act=act+X.all_states[pH][state]['pop']*self.kcat_microstates[state].get()
                    states=states+1
            act_prof1[pH]=act
        if states>0:
            max_act=max(act_prof1.values())
            if max_act==0.0:
                fact=1.0
            else:
                fact=1.0/max_act
            for ph in act_prof1.keys():
                act_prof1[ph] = fact*act_prof1[ph]
        else:
            act_prof1={}
        #
        # Secondary activity data
        #
        act_prof2={}
        states=0
        for pH in self.phvals:
            act=0.0
            for state in self.states:
                if self.act_state2[state].get()==1:
                    act=act+X.all_states[pH][state]['pop']
                    states=states+1
            act_prof2[pH]=act
        if states>0:
            max_act=max(act_prof2.values())
            if max_act==0.0:
                fact=1.0
            else:
                fact=1.0/max_act
            for ph in act_prof2.keys():
                act_prof2[ph] = fact*act_prof2[ph]
        else:
            act_prof2={}
        #
        # ftir data
        #
        if getattr(self,'FTIR_win',None):
            ftir = self.FTIR_win.get_sum_of_tits()
        else:
            ftir = 0
        
        return act_prof1,curve,ftir,act_prof2
    #
    #--------
    #

    def activity_diff(self, just_activity_fit = 0):
        act_prof_and_tc = self.get_act_prof_and_tc_and_ftir()
        diff = 0.0
        for ph in self.activity_data.keys():
            diff = diff + abs(act_prof_and_tc[0][ph]-self.activity_data[ph])

        if not just_activity_fit == 0: 
            if self.include_tcs_in_fit == 1:
                for group in self.exp_data.curves.keys():
                    for ph in self.exp_data.curves[group].keys():
                        diff = diff + abs(act_prof_and_tc[1].curves[group][ph]-self.exp_data.curves[group][ph])

        return diff
    #
    # ------
    #
    def evaluate_fit_ph_activity(self,make_win=1):
        import tkMessageBox

        if not getattr(self,'activity_data',None):
            tkMessageBox.showwarning('No pH-activity profile loaded',
                                     'Load pH-activity profile first')
            return

        no_active = 0
        for state in self.states:
            if self.act_state[state].get()==1:
                no_active = no_active+1

        if no_active == 0:
            tkMessageBox.showwarning('No active states selected',
                                     'Select active micro states first')
            return


        sum_unsigned_error = self.activity_diff(just_activity_fit=1)
        no_points = len(self.activity_data)
        mean_unsigned_error = sum_unsigned_error/no_points

        if make_win==1:
            text = """
            Sum of unsigned errors:    %8.2f
            Number of points:          %8d
            Mean unsigned error:       %8.4f""" %(sum_unsigned_error,no_points,mean_unsigned_error)

            tkMessageBox.showinfo('Evaluation of pH activity curve fit',
                                  text)

        return mean_unsigned_error

    #
    # -----
    #
    
    def stop_fit_activity(self):
        self.keep_running_activity = 0
        return

    #
    # -----
    #

    def cancel_copy_group_micro(self,event=None):
        """Cancel copy group to Ekin"""
        self.pick_group_micro.destroy()
        return
