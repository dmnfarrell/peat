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

from Tkinter import *
import Pmw

class make_canvas:

    def __init__(self,title,xmax,ymax,window,x_padding=10,y_padding=10):
        """Make a canvas"""
        self.sizex = min(max(xmax,400),1500)
	self.sizey = min(max(ymax,300),800)
	self.canvas_edge_x=self.sizex
        self.canvas_edge_y=self.sizey
        #
        # Calculate the x spacing - we base this on 100 frames - should be made dynamic
        # The X spacing cannot be lower than 1
        #
        self.x_padding=x_padding
        self.y_padding=y_padding
        
        #
        self.canv = Pmw.ScrolledCanvas(window,
                                       borderframe = 1,
                                       labelpos = 'n',
                                       label_text = title,
                                       usehullsize = 1,
                                       hull_width = self.canvas_edge_x,
                                       hull_height = self.canvas_edge_y,
                                       hscrollmode='dynamic',
                                       vscrollmode='dynamic'
                                       )
        self.canv.interior().configure(bg='white')
        self.set_scale()
        return

    #
    # ----
    #

    def get_x_spacing(self,steps):
        """Get the spacing if plotting steps items"""
        self.x_spacing=max(1,int(float(self.canvas_edge_x-self.x_padding)/float(steps)))
        return self.x_spacing

    def get_y_spacing(self,steps):
        """Get the spacing if plotting steps items"""
        self.y_spacing=max(1,int(float(self.canvas_edge_y-self.y_padding)/float(steps)))
        return self.y_spacing

    def get_x(self,x):
        return x*self.x_spacing+self.x_padding

    def get_y(self,y):
        """Get the y coordinate"""
        return self.sizey-self.y_padding-y*self.y_spacing
    
    def plot_value(self,x,y,value,colour=None):
        x_coord=self.get_x(x)
        y_coord=self.get_y(y)
        if not colour:
            colour=self.scale_colour(value)
        id=self.canv.create_rectangle(x_coord,y_coord,x_coord+self.x_spacing,y_coord+self.y_spacing,fill=colour,tags='box')
        return

    def set_x_label(self,label,x):
        """Put a label underneath the x-axis at x=x"""
        x=self.get_x(x+0.5)
        y=self.get_y(-1)
        text=''
        for char in label:
            text=text+char+'\n'
        id=self.canv.create_text(x,y,text=text,anchor='nw',tags='label')
        return
    
    def set_y_label(self,label,y):
        """Put a label underneath the x-axis at x=x"""
        x=self.get_x(-3)
        y=self.get_y(y-0.3)
        id=self.canv.create_text(x,y,text=label,anchor='nw',tags='label')
        return

    def set_scale(self,col_min=[0,0,255],col_mid=[255,255,255],col_max=[255,0,0],val_min=-1.0,val_mid=0.0,val_max=1.0):
        import numpy
        self.col_min=numpy.array(col_min)
        self.col_mid=numpy.array(col_mid)
        self.col_max=numpy.array(col_max)
        self.val_min=val_min
        self.val_mid=val_mid
        self.val_max=val_max
        return

    def scale_colour(self,value):
        """Scale a colour on our ramp"""
        if value is None:
            return 'black'
        import numpy
        value=min(self.val_max,value)
        value=max(self.val_min,value)
        if value<self.val_mid:
            colour=(value-self.val_min)/(self.val_mid-self.val_min)*(self.col_mid-self.col_min)+self.col_min
        else:
            value=abs(value)
            colour=(value-self.val_mid)/(self.val_max-self.val_mid)*(self.col_max-self.col_mid)+self.col_mid
        return self.colour_2_hex(colour)
         
    #
    # ----
    #

    def colour_2_hex(self,color):
        """Converts a color RGB tuple to a hex color"""
        r=hex(color[0]).replace('L','')
        g=hex(color[1]).replace('L','')
        b=hex(color[2]).replace('L','')
        return '#%s%s%s' %(r.split('x')[1].zfill(2),g.split('x')[1].zfill(2),b.split('x')[1].zfill(2))


class pKa_effect(Frame):

    def __init__(self,params):
        #
        # Initialise
        #
        self.params=params
        import pKaTool.pKaIO as pKaIO
        self.pdbfile=self.params['pdb']
        self.MCsteps=int(self.params['MCsteps'])
        self.pHstep=float(self.params['pHstep'])
        self.pHstart=float(self.params['pHstart'])
        self.pHend=float(self.params['pHend'])
        IO=pKaIO.pKaIO(self.pdbfile)
        #
        # Define the tmpdir for all our calcs
        #
        import os
        self.topdir=os.path.split(self.pdbfile)[0]
        self.topdir=os.path.join(self.topdir,'%s_autonomy' %os.path.split(self.pdbfile)[1])
        if not os.path.isdir(self.topdir):
            os.mkdir(self.topdir)
        #
        # Did we find a completed pKa calculation for the PDB file
        #
        if IO.calculation_completed:
            #
            # All files found
            #
            print 'pKa calculation files found'
        else:
            print
            print 'I could not find a completed pKa calculation for the specified PDB file'
            print 'Please complete the pKa calculation first'
            print
            raise Exception()
        #
        # Get the wild type titration curves calculated with WHAT IF
        #
        import pKaTool.pKaIO
        IO=pKaTool.pKaIO.pKaIO(self.pdbfile)
        self.wt_titcurv=IO.read_titration_curve()
        self.wt_pKas=IO.readpka()
        #
        # Recalculate titration curves with the CPP algorithm
        #
        import os
        dirname=os.path.join(self.topdir,'wt_recalc')
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        filename=os.path.join(dirname,'WT.DAT')
        data=None
        if os.path.isfile(filename):
            try:
                fd=open(filename)
                import cPickle
                data=cPickle.load(fd)
                fd.close()
                self.wtpkas=data['pkas']
                self.wt_titcurv=data['titcurv']
                print 'Loaded wild type recalculated titration curves'
                print 'Loaded %d titration curves' %(len(self.wt_titcurv.keys()))
            except:
                data=None
        #
        # If we don't have the data calculate it
        #
        if data is None:
            import pKa_MC
            print 'Recalculating wild type titration curves'
            pKa_params={'pHstart':self.pHstart,'pHstop':self.pHend,'pHstep':self.pHstep,'pKMCsteps':self.MCsteps,'verbose':1}
            self.pKaCALC=pKa_MC.pKa_calculation_class(self.pdbfile,pKa_info=None,params=pKa_params,parent=self)
            self.pKaCALC.set_MC_CPP()
            self.pKaCALC.set_reporter_groups(self.wt_pKas.keys())
            self.wtpkas,self.wt_titcurv=self.pKaCALC.calc_wt_pkas()
            #
            # Save the data
            #
            fd=open(filename,'w')
            import cPickle
            data=cPickle.dump({'pkas':self.wtpkas,'titcurv':self.wt_titcurv},fd)
            fd.close()
            print 'done'
            print
            print 'Calculated %d titration curves' %(len(self.wt_titcurv.keys()))
        #
        # Get the in_system energies
        # 
        self.in_data=self.calculate_insystem_dpKas()
        #
        # Get the ex_data
        #
        self.ex_data=self.calculate_exsystem_dpKas()
        return
    #
    # ----------------------
    #
            
    def load_groups(self):
        #
        # Get all the titratable groups in the pdb file
        #
        import Protool
        P=Protool.structureIO()
        P.readpdb(self.params['pdb'])
        self.groups=P.get_titratable_groups()
        return

    #
    # --------------------
    #

    def select_interface(self):
        """
        # Little interface to select the group(s) we're interested in
        """
        Frame.__init__(self)
        self.master.title('Analysing pKa value decomposition')
        rowspan=8
        Label(self.master,text='Select groups to see decomposition analysis').grid(row=0,column=0)
        row=1
        yscrollbar=Scrollbar(self.master,orient='vertical',width=10)
        yscrollbar.grid(row=row,column=2,rowspan=rowspan,sticky='nws')
        #
        height=10
        self.group_disp=Listbox(self.master,bg='white',
                            fg='black',
                            height=height,width=25,
                            yscrollcommand= yscrollbar.set,
                            selectmode=EXTENDED)
        self.group_disp.grid(row=row,column=0,columnspan=2,rowspan=rowspan,sticky='news')
        yscrollbar.config(command=self.group_disp.yview)
        self.group_disp.delete(0, END)
        for group in self.groups:
            self.group_disp.insert(END,group)
        #
        # Do the two canvases
        #
        sizex=600
        sizey=500
        self.insystem_canv=make_canvas('In_system',xmax=sizex,ymax=sizey,window=self.master,x_padding=125,y_padding=125)
        self.insystem_canv.canv.grid(row=0,rowspan=10,column=10,sticky='news')
        self.exsystem_canv=make_canvas('Ex_system',xmax=sizex,ymax=sizey,window=self.master,x_padding=125,y_padding=125)
        self.exsystem_canv.canv.grid(row=0,rowspan=10,column=21,sticky='news')
        self.diff_canv=make_canvas('ABS(In_system / Ex_system)',xmax=sizex,ymax=sizey,window=self.master,x_padding=125,y_padding=125)
        self.diff_canv.canv.grid(row=11,rowspan=10,column=10,sticky='news')
        #
        # Plot the two first matrices
        #
        self.plot_matrix(self.in_data,self.insystem_canv,tagname='insystem')
        self.plot_matrix(self.ex_data,self.exsystem_canv,tagname='exsystem')
        #
        # Construct the diff matrix
        #
        diff={}
        for group1 in self.groups:
            diff[group1]={}
            for group2 in self.groups:
                if group1==group2:
                    diff[group1][group2]=None
                elif self.ex_data[group1][group2]!=0.0:
                    ratio=self.in_data[group1][group2]/self.ex_data[group1][group2]
                    if ratio<0.0:
                        print '%10s-%10s in: %5.3f, ex: %5.3f, ratio: %5.3f' %(group1,group2,
                                                                           self.in_data[group1][group2],
                                                                           self.ex_data[group1][group2],
                                                                           ratio)
                                                                           
                        diff[group1][group2]='yellow'
                    else:
                        diff[group1][group2]=ratio
                elif self.in_data[group1][group2]!=0.0:
                    diff[group1][group2]='green'
                    print '%10s-%10s in: %5.3f, ex: %5.3f' %(group1,group2,
                                                                           self.in_data[group1][group2],
                                                                           self.ex_data[group1][group2])
                else:
                    diff[group1][group2]='grey'
        #
        # Plot it
        #
        print 'General legend'
        print '---------------'
        print '-1.0: red, 0.0: white, 1.0: blue'
        print
        print 'Additional legend for ratio plot'
        print 'Yellow: Ratio < 0.0'
        print 'White: Ratio = 1.0'
        print 'Red: ratio 0.0 -> 1.0 (white) -> 2.0 (blue)'
        print 'Grey: no dpKa'
        print 'Green: in_system value is 0.0, but ex_system value is not! - very interesting!!'
        print

        self.diff_canv.set_scale(val_mid=1.0,val_min=0.0,val_max=2.0)
        self.plot_matrix(diff,self.diff_canv,tagname='diff',diffmode=1)
        #
        # Write the CSV files
        #
        import os
        filename=os.path.join(self.topdir,'autonomy_result.csv')
        fd=open(filename,'w')
        fd.write('Removed group, Monitored group, dpKa_in, dpKa_ex\n')
        for r_group in self.groups:
            for m_group in self.groups:
                if r_group==m_group:
                    continue
                fd.write('%10s, %10s, %5.3f, %5.3f\n' %(r_group,m_group,self.in_data[r_group][m_group],self.ex_data[r_group][m_group]))
        fd.close()
        #
        # All Done
        #
        return

    #
    # ----
    #

    def plot_matrix(self,matrix,canvas,tagname,diffmode=None):
        """Plot the data in the matrix on the canvas"""
        groups=matrix.keys()
        groups.sort()
        xstep=canvas.get_x_spacing(len(groups))
        ystep=canvas.get_y_spacing(len(groups))
        #
        # Set the axis labels
        #
        x_end=len(groups)
        y_end=x_end
        canvas.set_x_label('Removed group',x_end)
        canvas.set_y_label('Monitored group',y_end)
        for removed_group in groups:
            canvas.set_x_label(removed_group,groups.index(removed_group))
            canvas.set_y_label(removed_group,groups.index(removed_group))
            for group2 in groups:
                x=groups.index(removed_group)
                y=groups.index(group2)

                if not matrix[removed_group].has_key(group2):
                    value=None
                else:
                    value=matrix[removed_group][group2]
                plotted=None
                if diffmode:
                    if type(value) is type(' '):
                        canvas.plot_value(x,y,value,colour=value)
                        plotted=1
                if not plotted:
                    canvas.plot_value(x,y,value)
        return

    #
    # --------------------
    #

    def calculate_insystem_dpKas(self):
        """Calculate the dpKas for all residues"""
        import pKarun.pKa_utility_functions as pKa_utility_functions
        import os
        groups=self.wt_pKas.keys()
        groups.sort()
        self.in_data={}
        count=1
        for group in groups:
            print 'Group: %s (%3d of %3d)' %(group,count,len(groups))
            count=count+1
            dirname=os.path.join(self.topdir,'in_system')
            if not os.path.isdir(dirname):
                os.mkdir(dirname)
            filename=os.path.join(dirname,group+'.DAT')
            data=None
            if os.path.isfile(filename):
                try:
                    fd=open(filename)
                    import cPickle
                    data=cPickle.load(fd)
                    fd.close()
                    self.in_data[group]=data
                    print 'Loaded in_system data for %s' %group
                except:
                    data=None
            #
            # If we don't have the data calculate it
            #
            if data is None:
                print 'Calculating in_system dpKas when removing: %s' %group
                import pKa_MC
                pKa_params={'pHstart':self.pHstart,'pHstop':self.pHend,'pHstep':self.pHstep,'pKMCsteps':self.MCsteps,'verbose':1}
                self.pKaCALC=pKa_MC.pKa_calculation_class(self.pdbfile,pKa_info=None,params=pKa_params,parent=self)
                self.pKaCALC.set_MC_CPP()
                reporter_groups=groups[:]
                reporter_groups.remove(group)
                self.pKaCALC.set_reporter_groups(reporter_groups)
                self.pKaCALC.remove_interactions_with(group)
                pkas,titcurvs=self.pKaCALC.calc_wt_pkas()
                self.in_data[group]=titcurvs
                print 'I got %d titration curves' %(len(titcurvs.keys()))
                #
                # Save the data
                #
                fd=open(filename,'w')
                import cPickle
                data=cPickle.dump(titcurvs,fd)
                fd.close()
                print 'done'
                print
            print
        #
        # Integrate differences from wild type titration curve
        #
        data={}
        rgroups=self.in_data.keys()
        rgroups.sort()
        for removed_group in rgroups:
            data[removed_group]={}
            tgroups=self.in_data[removed_group].keys()
            tgroups.sort()
            for tgroup in tgroups:
                if removed_group==tgroup:
                    continue
                pH_values=self.in_data[removed_group][tgroup].keys()
                pH_values.sort()
                pH_values.remove('pKa')
                pH_step=pH_values[1]-pH_values[0]
                diff=0.0
                for pH in pH_values:
                    diff=diff+(self.in_data[removed_group][tgroup][pH]-self.wt_titcurv[tgroup][pH])*pH_step
                if abs(diff)<0.1:
                    diff=0.0
                data[removed_group][tgroup]=diff
        return data

    #
    # ----
    #

    def do_ex_system(self):
        """
        # Calculate the effects of titratable groups "ex_system"
        """
        import os
        filename=os.path.join(self.topdir,'ex_system.DAT')
        data=None
        if os.path.isfile(filename):
            try:
                print 'Loading ex_system data'
                fd=open(filename)
                import cPickle
                data=cPickle.load(fd)
                fd.close()
                return data
            except:
                data=None
        #
        # If we don't have the data calculate it
        #
        if data is None:
            import pKaTool.pKaIO
            define_pka=pKaTool.pKaIO.pKaIO(self.pdbfile) 
            #
            # pKa calc found, let;s continue
            #
            pKa_values=define_pka.readpka()
            #
            # read the WHAT IF pKa files
            #
            desolv=define_pka.read_desolv() 
            backgr=define_pka.read_backgr()
            matrix=define_pka.read_matrix() 
            #
            # Get all groups
            #
            titgrps= matrix.keys()          
            titgrps.sort()#sorts titgroups
            #
            # Start calculating pKa values for each pair
            #
            pKa={}
            titcurv={}
            counter=0
            for group1 in titgrps:
                print 'Ex_system for group %s (%3d of %3d)' %(group1,counter,len(titgrps))
                counter=counter+1
                for group2 in titgrps:
                    #
                    # Do we need to insert an empty dictionary in pKa[group1]?
                    #
                    if not pKa.has_key(group1):
                        pKa[group1]={}
                        titcurv[group1]={}
                    #
                    # Do not calculate if group1 is equal to group2
                    #
                    if group1==group2:
                        continue
                    #
                    # Do not calculate if we already have the result
                    #
                    if pKa.has_key(group2):
                        if pKa[group2].has_key(group1):
                            pKa[group1][group2]=pKa[group2][group1]
                            titcurv[group1][group2]=titcurv[group2][group1]
                            continue
                    #
                    # Do the pKa calculation
                    #
                    import pKaTool.pKa_calc
                    define_MC=pKaTool.pKa_calc.Boltzmann()
                    define_MC.desolv={group1:desolv[group1],group2:desolv[group2]}
                    define_MC.backgr={group1:backgr[group1],group2:backgr[group2]}
                    define_MC.matrix={group1:matrix[group1],group2:matrix[group2]}
                    pKa[group1][group2]=define_MC.calc_pKas(0,0.1,-20.0,30.0)
                    titcurv[group1][group2]=define_MC.prot_states
            #
            # Save results
            #
            import pickle
            fd=open(filename,'w')
            pickle.dump(titcurv,fd)
            fd.close()
            print            
            print 'All done'
            print
        return titcurv

    #
    # -----
    #
    
    def calculate_exsystem_dpKas(self):
        """Read the ex_system file and calculate the dpKas"""
        #
        # Calculate the ex_system values
        #
        self.ex_tits=self.do_ex_system()
        #
        # Load the file
        #
        # We have the titration curves, now we just calculate differences to the model pK titration curve
        #
        data={}
        groups=self.wt_pKas.keys()
        groups.sort()
        for removed_group in groups:
            data[removed_group]={}
            for tgroup in groups:
                if tgroup==removed_group:
                    continue
                tit_in_presence=self.ex_tits[removed_group][tgroup][tgroup]
                del tit_in_presence['pKa']
                pH_values=tit_in_presence.keys()
                pH_values.sort()
                tit_expresence=self.get_intpka_titcurv(tgroup,pH_values=pH_values) 
                pH_step=pH_values[1]-pH_values[0]
                diff=0.0
                for pH in pH_values:
                    diff=diff+(tit_expresence[pH]-tit_in_presence[pH])*pH_step
                if abs(diff)<0.1:
                    diff=0.0
                data[removed_group][tgroup]=diff
        return data

    #
    # -----
    #

    def get_intpka_titcurv(self,tgroup,pH_values):
        """Return the titration curve of tgroup using pH_values and the model pKa value of the group"""
        modelpK=self.wt_pKas[tgroup]['modelpK']
        intpKa=modelpK+self.wt_pKas[tgroup]['desolv']+self.wt_pKas[tgroup]['backgr']
        acidbase=self.wt_pKas[tgroup]['acidbase'] #1 for bases, -1 for acids
        tc={}
        pH_values.sort()
        for pH in pH_values:
            if acidbase==-1:
                tc[pH]=-(1.0-(1.0/(1.0+10.0**(pH-intpKa))))
            else:
                tc[pH]=1.0/(1.0+10.0**(pH-intpKa))
        return tc
        
    #
    # --------------------
    #


    def analyse_results(self,dpKas,selected):
        #
        # Print some simple stats for now
        # Print delta pKa values
        #
        sel_list=[]
        for group in selected.keys():
            if selected[group]:
                sel_list.append(group)
        #
        # --
        #
        removed_groups=dpKas.keys()
        removed_groups.sort()
        data={}
        #
        # This way of adding effects of removing groups only works well for
        # a single observed group
        #
        for group in removed_groups:
            #
            # Find the dpKa for our group
            #
            text='%13s' %(self.groups[group])
            text='Effect of: %13s: ' %(self.groups[group])
            dpKa_sum=0.0
            for obs_group in sel_list:
                #print obs_group, group, self.groups[obs_group]
                #k=dpKas[group].keys()
                #k.sort()
                #print k
                #text=text+'%13s dpKa %6.3f' %(self.groups[obs_group],dpKas[group][self.groups[obs_group]])
                text=text+'%6.3f' %(dpKas[group][self.groups[obs_group]])
                dpKa_sum=dpKa_sum+dpKas[group][self.groups[obs_group]]
                #
                # If we observe more than one group, then we get the net shift, not the total shift
                # due to each group
                #
            print text
            #
            # Get the uniqueid for this group
            #
            data[self.groups[group]]=dpKa_sum
        #
        # Prepare the pie chart
        #
        real_dpKas={}
        for obs_group in sel_list:
            real_dpKas[self.groups[obs_group]]=self.wt_pKas[self.groups[obs_group]]['delec']
        #
        # Go
        #
        data,obs_group_text,title2=self.prepare_pie_chart(data,sel_list,real_dpKas)
        #
        # plot
        #
        import dislin_driver, os
        os.system('/usr/bin/eog %s' %(dislin_driver.pie_chart(data,'pie.tif',obs_group_text,title2)))
        #
        # Make a nice 2D-colour plot for all groups
        #
        matrix={}
        #
        # Construct the matrix. We have a value for all combinations
        #
        for group in removed_groups:
            if not matrix.has_key(group):
                matrix[group]={}
            for obs_group in removed_groups:
                matrix[group][obs_group]=abs(dpKas[group][self.groups[obs_group]])
        #
        # Plot it
        #
        import dislin_driver
        #os.system('/usr/bin/eog %s' %(dislin_driver.colour_2D(matrix,'in-system','HEWL','Removed group','Observed group','abs(dpKa)')))
        return
    
    #
    # -----------------------
    #

    def prepare_pie_chart(self,dict,sel_list,real_dpKas):
        #
        # Reformat the dict - we only display the max most important contributions
        #
        max=6 # Number of groups to display individually in pie chart
        vals=dict.values()
        abs_vals=[]
        for val in vals:
            abs_vals.append(abs(val))
        abs_vals.sort()
        #
        # Find 'max' biggest absolute values
        #
        newdict={}
        loop=len(abs_vals)
        print 'Number of dpKa values in total: %d' %loop
        #print abs_vals
        if loop>max:
            loop=max
        for x in range(loop):
            val=abs_vals[-1-x] # Get the biggest values (those are last in the list)
            name=None
            for key in dict.keys():
                if dict[key]==val or dict[key]==-val:
                    name=key
                    break
            if not name:
                raise 'Cannot find group'
            newdict[name]=dict[key]
            del dict[key]
        #
        # How much does the rest contribute?
        #
        # Get sum of most important groups
        #
        newsum=0.0
        for key in newdict.keys():
            #print key,newdict[key]
            newsum=newsum+newdict[key]
        #print 'Most important residues',newdict.keys()
        #
        # Get sum of all the other groups
        #
        oldsum=0.0
        residues=dict.keys()
        residues.sort()
        #print 'All residues',residues
        for key in dict.keys():
            oldsum=oldsum+dict[key]
        rest=oldsum
        omitted=len(dict.keys())
        newdict['Remaining: %3d groups' %(omitted)]=rest
        #
        # Prepare info on obs_groups
        #
        obs_group_text=''
        for o_group in sel_list:
            obs_group_text=obs_group_text+self.groups[o_group]+' '
        #
        # Calculate difference between sum of decomposed dpKas and real dpKa
        #
        sum_diff=0.0
        for o_group in sel_list:
            sum_diff=sum_diff+abs(real_dpKas[self.groups[o_group]])
        text='Sum of real dpKas: %6.3f, ' %sum_diff
        text=text+'\nSum of decomposed dpKas: %6.3f, ' %(oldsum+newsum)
        text=text+'\nDiff.: %6.3f' %((oldsum+newsum)-sum_diff)
        print text
        #
        # Done
        #
        return newdict,obs_group_text,text
#
# ------------------------------
#


def parse_arguments(arg_list):
    #
    # First we define the defaults for the program, then
    # we parse the arguments
    # -------------
    #
    # Define defaults
    #
    # Types for parameter passing:
    # text: text. res_num_list: List of residu numbers + associated number for each
    # number: number
    # res_list: List of residue numbers
    # T/F: True/False: If you give the argument then set to true, otherwise false
    # 
    defaults={'pdb':[None,'text'],'pKas':[None,'res_num_list'],'pHstep':[0.01,'number'],
    'pHstart':[-20.0,'number'],'pHend':[20.0,'number'],'MCsteps':[200000,'number']}
    #
    # Parse the arguments
    #
    import console_interface
    parms=console_interface.parse_arguments(arg_list,defaults)
    # Reformat
    realparms={}
    for key in parms.keys():
        realparms[key]=parms[key][0]
    if not realparms['pdb']:
        print
        print 'You have to give a pdb file name with -pdb <filename>'
        print
        raise Exception()
    return realparms

#
# ------------------------------
#

def main():
    #
    # Parse arguments
    #
    import sys
    params=parse_arguments(sys.argv)
    #
    # Instantiate our little class
    #
    I=pKa_effect(params)
    #
    # Load all the groups
    #
    I.load_groups()
    #
    # Select group of interest
    #
    I.select_interface()
    I.mainloop()
    #
    # Calculate the results
    #
    #if choice==1:
    #    dpKas=I.calculate_dpkas(I.in_system)
    #    I.analyse_results(dpKas,selected)
    #elif choice==2:
    #    dpKas=I.calculate_dpkas_fromdict(I.ex_system)
    #    I.analyse_results(dpKas,selected)
    #
    # Done
    #
    return

if __name__=='__main__':
    main()
