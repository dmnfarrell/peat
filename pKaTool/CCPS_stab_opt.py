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
#
from Tkinter import *
from numpy import *
from numpy.linalg import * 
inverse=inv
Float=float

from titration_class import *
    
import math

class Optimisation_Analysis:

    def stab_and_CCPS_pop(self):
        """Driver for optimising stability and CCPS population"""
        #
        # Open the optimisation window
        #
        self.optwin=Toplevel()
        parent_geom=self.window.winfo_geometry()
        parent_geom=parent_geom.split('+')[1:]
        # Position optimisation window at same position as parent (titration curve window
        self.optwin.geometry('+%d+%d' %(int(parent_geom[0])+100,int(parent_geom[1])+200))
        self.optwin.title('System optimisation and analysis')
        #
        # Add scoring function specification interface
        #
        Label(self.optwin,text='Name').grid(row=1,column=0)
        Label(self.optwin,text='Optimise').grid(row=1,column=1)
        Label(self.optwin,text='Opt type').grid(row=1,column=2)
        Label(self.optwin,text='Cutoff value').grid(row=1,column=3)
        Label(self.optwin,text='Weight').grid(row=1,column=4)
        Label(self.optwin,text='pH range').grid(row=1,column=5)
        Label(self.optwin,text='Use profile').grid(row=1,column=6)
        Label(self.optwin,text='Optimise sum').grid(row=1,column=7)
        #
        self.properties={'CCPS population':{},
                         'CCPS2 population':{},
                         'CCPS maxrange':{},
                         'Stability':{},
                         'individual intpKa shift':{},
                         'individual intene':{}
                         }
        props=self.properties.keys()
        props.sort()
        row=1
        for prop in props:
            row=row+1
            Label(self.optwin,text=prop).grid(row=row,column=0)
            #
            self.properties[prop]['include']={'var':IntVar()}
            self.properties[prop]['include']['widget']=Checkbutton(self.optwin,variable=self.properties[prop]['include']['var'],onvalue=1,offvalue=0)
            self.properties[prop]['include']['widget'].grid(row=row,column=1)
            self.properties[prop]['include']['var'].set(1)
            #
            self.properties[prop]['opttype']={'var':StringVar()}
            self.properties[prop]['opttype']['widget']=Menubutton(self.optwin,textvariable=self.properties[prop]['opttype']['var'],relief='raised')
            self.properties[prop]['opttype']['menu']=Menu(self.properties[prop]['opttype']['widget'],tearoff=0)
            (self.properties[prop]['opttype']['widget'])['menu']=self.properties[prop]['opttype']['menu']
            #
            # Set defaults
            #
            if prop=='CCPS population':
                value=90.0
                defopttype='maximise'
                weight=1.0
            elif prop=='CCPS2 population':
                value=10.0
                defopttype='maximise'
                weight=1.0     
            elif prop=='Stability':
                value=0.0
                defopttype='maximise'
                weight=1.0
            elif prop=='individual intpKa shift':
                value=2.0
                defopttype='keep below'
                weight=10000.0
            elif prop=='individual intene':
                value=4.6
                defopttype='keep below'
                weight=10000.0
            elif prop=='CCPS maxrange':
                value=0.20
                defopttype='maximise'
                weight=1.0
            else:
                raise Exception
            #
            # Add all opttypes
            #
            for opttype in ['maximise','minimise','keep at','keep above','keep below']:
                self.properties[prop]['opttype']['menu'].add_radiobutton(label=opttype,
                                                                         variable=self.properties[prop]['opttype']['var'],
                                                                         value=opttype,
                                                                         indicatoron=1)
            self.properties[prop]['opttype']['widget'].grid(row=row,column=2)
            self.properties[prop]['opttype']['var'].set(defopttype)
            #
            self.properties[prop]['cutoff']={'var':DoubleVar()}
 
            self.properties[prop]['cutoff']['var'].set(value)
            self.properties[prop]['cutoff']['widget']=Entry(self.optwin,textvariable=self.properties[prop]['cutoff']['var'],width=10)
            self.properties[prop]['cutoff']['widget'].grid(row=row,column=3)
            #
            self.properties[prop]['weight']={'var':DoubleVar()}
            self.properties[prop]['weight']['var'].set(weight)
            self.properties[prop]['weight']['widget']=Entry(self.optwin,textvariable=self.properties[prop]['weight']['var'],width=10)
            self.properties[prop]['weight']['widget'].grid(row=row,column=4)
            #
            self.properties[prop]['pH']={'var':StringVar()}
            self.properties[prop]['pH']['var'].set('6.5-7.5')
            self.properties[prop]['pH']['widget']=Entry(self.optwin,textvariable=self.properties[prop]['pH']['var'],width=10)
            self.properties[prop]['pH']['widget'].grid(row=row,column=5)
            #
            # Add a profile selector
            #
            if prop in ['CCPS maxrange','CCPS population','CCPS2 population','Stability']:
                self.properties[prop]['profile']={'var':IntVar()}
                self.properties[prop]['profile']['var'].set(0)
                self.properties[prop]['profile']['widget']=Checkbutton(self.optwin,
                    variable=self.properties[prop]['profile']['var'],onvalue=1,offvalue=0)
                self.properties[prop]['profile']['widget'].grid(row=row,column=6)
            #
            # Selector for individual opt of protonation state populations
            #
            if prop in ['CCPS maxrange','CCPS population','CCPS2 population']:
                self.properties[prop]['sumopt']={'var':IntVar()}
                self.properties[prop]['sumopt']['var'].set(1)
                self.properties[prop]['sumopt']['widget']=Checkbutton(self.optwin,
                                                                       variable=self.properties[prop]['sumopt']['var'],onvalue=1,offvalue=0)
                self.properties[prop]['sumopt']['widget'].grid(row=row,column=7)
            
        #
        # Open the microscopic population states
        #
        self.micro_var.set(1)
        self.update_pkasystem_curves(0)
        #
        # Open the stability window
        #
        self.stab_button.select()
        self.stability_on_off()
        #
        # Continue with the optimisation window
        #
        row=row+1
        Label(self.optwin,text='CCPS chosen').grid(row=row,column=0)
        self.chosen_CCPS=StringVar()
        Label(self.optwin,textvariable=self.chosen_CCPS).grid(row=row,column=1)
        self.find_chosen_CCPS()
        Button(self.optwin,text='update',command=self.find_chosen_CCPS).grid(row=row,column=2)
        #
        # Add action buttons
        #
        row=row+1
        Button(self.optwin,text='Load primary pH-activity profile',command=self.pre_load_pH_activity_profile).grid(row=row,column=0)
        Button(self.optwin,text='Load second pH-activity profile',command=self.pre_load_2nd_pH_activity_profile).grid(row=row,column=1)
        Button(self.optwin,text='Load pH-stability profile',state=DISABLED).grid(row=row,column=2)
        #
        row=row+1
        Button(self.optwin,text='Randomize system',fg='pink',command=self.randomizesystem).grid(row=row,column=0)
        Button(self.optwin,text='Start optimisation',fg='green',command=self.optimise_system).grid(row=row,column=1)
        Button(self.optwin,text='Stop optimisation',fg='red',command=self.stopopt).grid(row=row,column=2)
        
        #
        # Status bar/window
        #
        row=row+1
        self.status=StringVar()
        Label(self.optwin,textvariable=self.status,bg='white').grid(row=row,column=0,columnspan=4)
        #
        # Adjust the resolution of sliders
        #
        for group in self.groups:
            self.groups[group].intpka_scale.configure(resolution=0.01)
            for group2 in self.groups[group].intenes.keys():
                self.groups[group].intenes_scales[group2].configure(resolution=0.01)
        #
        # Set pH step
        #
        self.pHstep.set(0.1)
        return

    #
    # ----
    #

    def pre_load_pH_activity_profile(self):
        """Handle the pH_activity profile load"""
        self.load_pH_activity_profile(parent=self.optwin)
        self.properties['CCPS population']['profile']['var'].set(1)
        return

    def pre_load_2nd_pH_activity_profile(self):
        """Handle the pH_activity profile load"""
        self.load_pH_activity_profile(parent=self.optwin)
        self.properties['CCPS2 population']['profile']['var'].set(1)
        return

    #
    # ----
    #

    def stopopt(self):
        """Stop the current fitting operation"""
        self.keep_running_opt=0
        return

    #
    # -----
    #

    def find_chosen_CCPS(self):
        """Sets the variable that displays the protonation states that are chosen for the CCPS"""
        self.chosen_CCPS.set(self.get_CCPSs())
        return
    
    #
    # -----
    #

    def randomizesystem(self):
        return

    #
    # -----
    #

    def optimise_system(self):
        """Optimise the system according to the values given by the user"""
        #
        # Did the bozo choose the CCPS?
        #
        self.find_chosen_CCPS()
        if self.get_CCPSs()==[]:
            import tkMessageBox
            tkMessageBox.showwarning('No CCPS selected',
                                     'You have to chose the CCPS from the microscopic states',
                                     parent=self.optwin)
            return
        #
        # Find the properties we will use
        #
        self.active_properties=[]
        props=self.properties.keys()
        for prop in props:
            if self.properties[prop]['include']['var'].get()==1:
                self.active_properties.append(prop)
        #
        print 'I am optimising these properties:',self.active_properties
 
        #
        # Get the variables
        #
        self.intpkas=[]
        self.intenes=[]
        for group in self.groups:
            self.intpkas.append(self.groups[group].intpka)
            for group2 in self.groups[group].intenes.keys():
                self.intenes.append(self.groups[group].intenes[group2])
        #
        # Perform all variations within the limits
        #
        self.intpka_limit=self.properties['individual intpKa shift']['cutoff']['var'].get()
        self.intene_limit=self.properties['individual intene']['cutoff']['var'].get()
        print 'Intpka limit: %5.3f' %self.intpka_limit
        print 'intene limit: %5.3f' %self.intene_limit
        #
        # Make copies of all the variables
        #
        self.vars=[]
        for group in self.groups:
            self.vars.append([self.groups[group].intpka,'intpka',group])
            for group2 in self.groups[group].intenes.keys():
                self.vars.append([self.groups[group].intenes[group2],'intene',group])
        #
        # Set damper
        #
        self.LM_damper = 0.1
        #
        # Update the counter
        #
        self.count=0
        self.keep_running_opt=1
        #
        # Init vars for the scoring function
        #
        self.CCPS_pH_vals=None
        self.CCPSmax_pH_vals=None
        self.Stability_pH_vals=None
        #
        # Start off the iterations
        #
        scores=self.get_system_score().values()
        old_diff=10000.0
        now_diff=0.0
        for val in self.get_system_score().values():
            now_diff=now_diff+val
        min_steps=100
        for x in range(1,10000):
            self.fit_LM_optimise()
            #
            # Update the counter, the scales and the titration curves
            #
            self.status.set("Optimisation running. Step %5d" %self.count)
            self.count=self.count+1
            self.window.update()
            self.optwin.update()

            #self.optwin.update_idletasks()
            self.update_scales_from_fit()
            self.titwin.update() # This statement is essential for propagating changes in variable values
            #
            #
            # recalculate the differences
            #
            now_diff=0.0
            for val in self.get_system_score().values():
                now_diff=now_diff+val
            #
            # Check convergence
            #
            if abs(now_diff-old_diff)<0.000005 and x>min_steps:
                print 'Converged',now_diff
                break
            else:
                old_diff=now_diff
            #
            # Should we stop?
            #
            if self.keep_running_opt==0:
                print 'User abort'
                break
        self.update_pkasystem_curves()
        return

    #
    # ----
    #

    def get_jacobian_optimise(self):
        """Get the Jacobian matrix and errors of the data points"""
        #
        # We have one score for each property included
        #
        no_data_points=len(self.active_properties)
        if 'CCPS2 population' in self.active_properties:
            no_data_points=no_data_points-1
        errors = resize(array(0,float),[no_data_points])        
        jacobian = resize(array(0,float),[no_data_points,len(self.vars)])
        #
        # Precalculate the variation of all parameters
        #
        
        variations=[]
        step = 1e-8
        for var in range(len(self.vars)):
            self.vars[var][0].set(self.vars[var][0].get()+step)
            variations.append(self.get_system_score())
            self.vars[var][0].set(self.vars[var][0].get()-step)
        #
        # construct jacobian
        #
        data_id=0
        #
        # Find the errors 
        #
        now=self.get_system_score()
        for prop in self.active_properties:
            if prop=='CCPS2 population':
                continue
            y=0 # All functions always have minimum at zero
            errors[data_id] =0-now[prop]
            #
            # Find the derivatives
            #
            diff=resize(array(0,float),[len(self.vars)])
            count=0
            for variation in variations:
                diff[count]=(now[prop]-variation[prop])/step
                count=count+1
            jacobian[data_id]=diff
            data_id=data_id+1
        return jacobian,errors

    #
    # ----
    #

    def calc_sq_diff_optimise(self,charges):
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

    def fit_LM_optimise(self):
        """Do Levenberg-Marquardt fitting"""
        J,E =self.get_jacobian_optimise()
        JT = transpose(J)
        JTE = dot(JT,E)
        JTJ = dot(JT,J)
        JTJd = JTJ + self.LM_damper*identity(shape(JTJ)[0])
        invJTJd = inverse(JTJd)
        q = -dot(JTE,invJTJd)
        #
        # Change the values
        #
        for var in range(len(self.vars)):
            #
            # Update the vars if they do not get outside the bounds
            #
            new_val=self.vars[var][0].get()+(q[var]/abs(q[var]))*min(abs(q[var]),0.2)
            group=self.vars[var][2]
            ok=None
            if self.vars[var][1]=='intpka':
                dintpka=abs(new_val-self.unfolded_groups[group].intpka.get())
                if dintpka<=self.intpka_limit and new_val>=0.0:
                    ok=1
                if self.properties['individual intpKa shift']['include']['var'].get()==0:
                    #
                    # If no restriction then the shift is always ok
                    #
                    ok=1
            elif self.vars[var][1]=='intene':
                if new_val<=self.intene_limit and new_val>=0.0:
                    ok=1
                if self.properties['individual intene']['include']['var'].get()==0:
                    #
                    # If no restriction then the shift is always ok
                    #
                    ok=1
            else:
                print 'Unknown variable'
                raise Exception
            #
            # If ok, then update vars
            #
            if ok:
                self.vars[var][0].set(new_val)
        return

    #
    # -----
    #
    
    def get_system_score(self):
        """Get value of the scoring function specified by the user"""

        self.X,pKa_values=self.calc_pKas(self.groups)
        curve=titration_curve(self.X.prot_states)
        self.stability=self.do_stab_curve(self.X)


        self.act_prof={}
        for pH in self.phvals:
            act=0.0
            for state in self.states:
                if self.act_state[state].get()==1:
                    act=act+self.X.all_states[pH][state]['pop']
            self.act_prof[pH]=act
        #
        # Secondary activity profile
        #
        self.act_prof2={}
        for pH in self.phvals:
            act=0.0
            for state in self.states:
                if self.act_state2[state].get()==1:
                    act=act+self.X.all_states[pH][state]['pop']
            self.act_prof2[pH]=act
        #
        # Score everything
        #
        scores={}
        CCPS2_score=None
        for prop in self.active_properties:
            thisprop=self.properties[prop]
            function=getattr(self,'_'+prop.replace(' ','_'))
            if prop=='CCPS2 population':
                CCPS2_score=function(thisprop)
            else:
                scores[prop]=function(thisprop)
        #
        # Score CCPSs under one
        #
        if CCPS2_score:
            scores['CCPS population']=max(CCPS2_score,scores['CCPS population'])
            print 'Doing max'
        #
        #
        props=scores.keys()
        props.sort()
        for prop in props:
            print '%10s score: %9.5f,' %(prop,scores[prop])
        print
        return scores

    #
    # -----
    #

    def _CCPS_population_base(self,variable,primary,spec):
        """Function that will give a score for how well one of the CCPS population meets the targets"""
        #
        # Loaded profile or entered values
        #
        if spec['profile']['var'].get()==0:
            #
            # Entered values
            #
            # Get the pH values where we want to impose the criterium
            #
            if not self.CCPS_pH_vals:
                self.CCPS_pH_vals=self.get_pH_range(spec['pH']['var'].get())
            print 'Monitoring activity at these pH values',self.CCPS_pH_vals
            #
            # pH values have already been identified
            #
            opttype=spec['opttype']['var'].get()
            cutoff=spec['cutoff']['var'].get()
            if opttype=='maximise':
                sum_act=0.0
                count=0
                for pH in self.CCPS_pH_vals:
                    if primary:
                        if self.act_prof[pH]<cutoff:
                            sum_act=sum_act+self.act_prof[pH]
                            count=count+1
                    else:
                        if self.act_prof2[pH]<cutoff:
                            sum_act=sum_act+self.act_prof2[pH]
                            count=count+1
            
                if count>1:
                    sum_act=sum_act/float(count)
                sum_act=max(cutoff/100.0-sum_act,0)
            else:
                raise Exception
        else:
            #
            # Loaded profile
            #
            sum_act=0.0
            #
            # Optimse profile, or each protonation state on its own
            #
            if spec['sumopt']['var'].get()==1:
                for pH in variable.keys():
                    if primary:
                        sum_act=sum_act+abs(variable[pH]-self.act_prof[pH])
                    else:
                        sum_act=sum_act+abs(variable[pH]-self.act_prof2[pH])
            else:
                #
                # Find the difference between each protonation state
                #
                diffs=[]
                for state in self.states:
                    thisdiff=0.0
                    if (primary and self.act_state[state].get()==1) or (not primary and self.act_state2[state].get()==1):
                        for pH in variable.keys():
                            thisdiff=thisdiff+abs(self.X.all_states[pH][state]['pop']-variable[pH])
                        diffs.append(thisdiff)
                sum_act=max(diffs)
        #
        # Scale the sum by the weight
        #
        sum_act=sum_act*spec['weight']['var'].get()
        return sum_act

    #
    # ----
    #

    def _CCPS_population(self,spec):
        """Get the score for the primary CCPS population"""
        return self._CCPS_population_base(self.activity_data,1,spec)

    def _CCPS2_population(self,spec):
        """Get the score for the seconary CCPS population"""
        return self._CCPS_population_base(self.secondary_activity_data,1,spec)


    #
    # ----
    #

    def _CCPS_maxrange(self,spec):
        """Score for having the maximum of the CCPS in a range"""
        #
        # Get the pH values where we want to impose the criterium
        #
        if not self.CCPSmax_pH_vals:
            self.CCPSmax_pH_vals=self.get_pH_range(spec['pH']['var'].get())
            print 'Monitoring activity at these pH values',self.CCPSmax_pH_vals
        #
        # pH values have already been identified
        #
        opttype=spec['opttype']['var'].get()
        cutoff=spec['cutoff']['var'].get()
        if opttype=='maximise':
            #
            # Find the maximum activity
            #
            maxact=0.0
            for pH in self.act_prof.keys():
                maxact=max(maxact,self.act_prof[pH])
            #
            # Report the fraction of points that are within 10%
            #
            count=0
            ok=0.0
            for pH in self.CCPSmax_pH_vals:
                if maxact>0.000001:
                    if abs((self.act_prof[pH]-maxact)/maxact)<=cutoff:
                        pass
                    else:
                        ok=ok+abs((self.act_prof[pH]-maxact)/maxact)
                        print pH,(self.act_prof[pH]-maxact)/maxact
                else:
                    ok=ok+10*(0.000001-maxact)
                count=count+1
            if count>0:
                ok=ok/float(count)
            #
        else:
            raise Exception
        ok=ok*spec['weight']['var'].get()
        return ok

    #
    # ----
    #

    def _Stability(self,spec):
        """Score the stability of the system"""
        if not self.Stability_pH_vals:
            self.Stability_pH_vals=self.get_pH_range(spec['pH']['var'].get())
        #
        # Score the stability
        #
        opttype=spec['opttype']['var'].get()
        cutoff=spec['cutoff']['var'].get()
        if opttype=='maximise':
            # Maximising stability means minimising the dG of folding
            sum_stab=0
            count=0
            for pH in self.Stability_pH_vals:
                if self.stability[pH]>cutoff:
                    sum_stab=sum_stab+self.stability[pH]-cutoff
                    count=count+1
            if count>1:
                sum_stab=sum_stab/float(count)
        else:
            raise Exeption
        import math
        sum_stab=sum_stab*spec['weight']['var'].get()
        #print 'Stability',math.exp(sum_stab)
        return sum_stab

    #
    # ----
    #

    def _individual_intpKa_shift(self,spec):
        """Score the total pKa shift"""
        return 0.0
        opttype=spec['opttype']['var'].get()
        cutoff=spec['cutoff']['var'].get()
        if opttype=='keep below':
            score=0.0
            for group in self.groups:
                dintpka=abs(self.groups[group].intpka.get()-self.unfolded_groups[group].intpka.get())
                if dintpka>cutoff:
                    score=score+math.pow(dintpka-cutoff,4)
        else:
            raise Exception
        return score*spec['weight']['var'].get()

    #
    # ----
    #

    def _individual_intene(self,spec):
        return 0.0
        """Keep all interaction energies below a certain kT limit"""
        import math
        opttype=spec['opttype']['var'].get()
        cutoff=spec['cutoff']['var'].get()
        if opttype=='keep below':
            score=0.0
            for group in self.groups:
                for num in self.groups[group].intenes:
                    intene=self.groups[group].intenes[num].get()
                    if intene>cutoff:
                        score=score+math.pow(intene-cutoff,4)
        else:
            raise Exception
        return score*spec['weight']['var'].get()
    #
    # ----
    #


    def get_pH_range(self,pH_vals):
        """Return a list of pH values that should be monitored in the individual scoring functions
        from the input given in the pH range field"""
        ranges=pH_vals.split(',')
        pH_vals=[]
        for sep in ranges:
            if sep.find('-')!=-1:
                start_end=sep.split('-')
                start=float(start_end[0])
                end=float(start_end[1])
                for pH in self.act_prof.keys():
                    if pH>=start and pH<=end:
                        pH_vals.append(pH)
            else:
                #
                # Simple number - just pick a value closer than 0.1
                #
                pH_set=float(sep)
                if self.act_prof.has_key(pH_set):
                    pH_vals.append(pH_set)
                else:
                    for pH in self.act_prof.keys():
                        if abs(pH-pHset)<0.1:
                            pH_vals.append(pH)
        pH_vals.sort()
        return pH_vals
