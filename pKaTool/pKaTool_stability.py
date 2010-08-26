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
from pKarun.pKarun_base import *
import pKa_system

class pKaTool_stability(pKarun):

    #
    # -----
    #

    def do_stab_curve(self):
        """Calculate the pH-dependence of stability"""
        #
        # Calculate the stability curve
        #
        #
        # The pKa values for the unfolded form are simply the model pKa values
        #
        import math, pKarun.pKa_utility_functions
        UF={}
        for group in self.selected_groups:
            UF[group]={}
            model=self.pkavals[group]['modelpK']
            phvals=self.titration_curves[group].keys()
            phvals.sort()
            for pH in phvals[:-1]:
                crg=1.0/(math.pow(10.0,float(pH)-model)+1)
                if pKarun.pKa_utility_functions.isacid(group):
                    crg=-1+crg
                UF[group][pH]=crg
        #
        # Integrate - do it in two ways
        #
        integral=0.0
        integral2=0.0
        intcurve=[]
        intcurve2=[]
        group0=self.pkavals.keys()[0]
        pHvalues=self.titration_curves[group0].keys()
        pHvalues.sort()
        pHvalues=pHvalues[:-1]
        dpH=abs(pHvalues[0]-pHvalues[1])
        min_val=99999
        max_val=-9999
        import math
        ln10=math.log(10)
        #
        # Convert to kJ/mol - disabled
        #
        #k=1.3806503E-23
        #T=298.15
        #Na=6.02214199E23
        #factor=k*T*Na/1000.0
        factor=1
        #
        summed_contributions_precalc={}
        for pH in pHvalues:
            intcurve.append(integral)
            intcurve2.append(integral2)
            if self.contrib_type.get()==1:
                if integral<min_val:
                    min_val=integral
                if integral>max_val:
                    max_val=integral
            integral2=0
            for group in self.selected_groups:
                integral=integral+ln10*dpH*(self.titration_curves[group][pH]-UF[group][pH])*factor
                for group2 in self.selected_groups:
                    if group!=group2:
                        integral2=integral2+(abs(self.titration_curves[group][pH])*abs(self.titration_curves[group2][pH])*self.matrix[group][group2][0])/2.0
                #
                # Sum the contributions from each group - this is a precalculation
                #
                if not summed_contributions_precalc.has_key(group):
                    summed_contributions_precalc[group]=0.0
                if self.contrib_type.get()==1:
                    summed_contributions_precalc[group]=summed_contributions_precalc[group]+ln10*dpH*(self.titration_curves[group][pH]-UF[group][pH])*factor
                else:
                    stab=0.0
                    for group2 in self.selected_groups:
                        if group!=group2:
                            stab=stab+(abs(self.titration_curves[group][pH])*abs(self.titration_curves[group2][pH])*self.matrix[group][group2][0])/2.0
                    summed_contributions_precalc[group]=max(summed_contributions_precalc[group],abs(stab))

            #
            # Update the min and max
            #
            if integral2<min_val:
                min_val=integral2
            if integral2>max_val:
                max_val=integral2
        #
        # From the pre-calculated contributions for each indiviual residue - determine the
        # the best cutoff for viewing the contributions (max 20 groups in display)
        #
        if self.show_grp_contribs.get()==1:
            num_groups=int(self.stab_cutoff.get())
            show_groups=[]
            contribs=[]
            for group in self.selected_groups:
                contribs.append([abs(summed_contributions_precalc[group]),group])
            contribs.sort()
            show_groups=[]
            max_val=0.0
            min_val=0.0
            sum=0.0
            for group in contribs[-min(num_groups,len(contribs)):]:
                show_groups.append(group[1])
                sum=sum+abs(group[0])
            max_val=sum/2.0
            min_val=-sum/2.0
        #
        # Draw the lines
        #
        lastpH=pHvalues[0]
        lastval=intcurve[0]
        count=1
        span=max_val-min_val
        #
        # Delete the lines from last time
        #
        for line in self.stab_lines.keys():
            self.stab_tc.delete(line)
            del self.stab_lines[line]
        #
        # Draw the y axis
        #
        canvas=self.stab_tc
        x_axis=self.get_x(pHvalues[0])-20
        y_axis=pKa_system.get_y_fromstab(min_val,span)
        endy=pKa_system.get_y_fromstab(max_val,span)
        self.stab_lines[self.stab_tc.create_line(x_axis,max([160,y_axis]),
                                                  x_axis,endy-10,fill='black',
                                                  width=self.linewidth)]=1
        self.stab_lines[self.stab_tc.create_text(x_axis,endy-35,text='dG of folding (kT)',
                                           fill='black',
                                           anchor='w')]=1
        #self.stab_lines[self.stab_tc.create_text(x_axis,endy-55,text='Net Electrostatic interactions (folded - unfolded) (kT)',
        #                                   fill='red',
        #                                   anchor='w')]=1
        #
        # Tick marks and tick labels
        #
        for tickval in range(int(min_val*100),int(max_val*100),int(max([(span*100.0)/5.0,1.0]))):
            y=pKa_system.get_y_fromstab(tickval/100.0,span)
            self.stab_lines[self.stab_tc.create_line(x_axis,
                                                y,x_axis-5,y,
                                               fill='black',width=self.linewidth)]=1
            self.stab_lines[self.stab_tc.create_text(x_axis-25,y,text='%5.1f' %(
                float(tickval)/100.0),fill='black')]=1
        #
        # Draw the stability line
        #
        summed_contributions={}
        label_position={}
        count=1
        #
        # Loop over pH values
        #
        for pH in pHvalues[1:]:
            lastx=self.get_x(lastpH)
            lasty=pKa_system.get_y_fromstab(lastval,span)
            val=intcurve[count]
            count=count+1
            x=self.get_x(pH)
            y=pKa_system.get_y_fromstab(val,span)
            self.stab_lines[self.stab_tc.create_line(lastx,lasty,float(x),float(y),
                                                  fill='black',
                                                  width=self.linewidth)]=1
            #
            # Outline the contribution of each group
            #
            if self.show_grp_contribs.get()==1:
                colour_count=0
                null_y=pKa_system.get_y_fromstab(0.0,span)
                starty_positive=null_y
                starty_negative=null_y
                ufgroups=self.selected_groups
                ufgroups.sort()
                for group in ufgroups:
                    #
                    # We only show the group if the summed contribution is higher than the cutoff
                    #
                    if group in show_groups:
                        #
                        # Make sure the dictionary is initialised
                        #
                        if not summed_contributions.has_key(group):
                            summed_contributions[group]=0.0
                            label_position[group]=None
                        #
                        # Get this contribution 
                        #
                        dx=abs(lastx-x)
                        if self.contrib_type.get()==1:
                            #
                            # Here we get the stability contribution from pKa shifts
                            #
                            endy=pKa_system.get_y_fromstab(ln10*dpH*(self.titration_curves[group][pH]-UF[group][pH])*factor,span)-null_y
                            summed_contributions[group]=summed_contributions[group]+endy
                        else:
                            #
                            # Otherwise the stability contribution from charge-charge interactions
                            #
                            stab=0.0
                            for group2 in self.selected_groups:
                                #
                                # Get the interaction between this group and the other group
                                #
                                if group!=group2:
                                    stab=stab+abs(self.titration_curves[group][pH])*abs(self.titration_curves[group2][pH])*self.matrix[group][group2][0]/2.0*factor
                            #
                            # Sum it
                            #
                            endy=pKa_system.get_y_fromstab(stab,span)-null_y
                            summed_contributions[group]=endy
                        #
                        # Draw the box
                        #
                        endy=summed_contributions[group]
                        if endy>0:
                            self.stab_lines[self.stab_tc.create_rectangle(x+1.5*dx,starty_positive,lastx+1.5*dx,endy+starty_positive,
                                                                          fill=self.colour_order[colour_count],
                                                                          outline='',
                                                                          stipple='gray50',
                                                                          width=self.linewidth)]=1
                            label_position[group]=(starty_positive*2+endy)/2.0
                            starty_positive=endy+starty_positive
                        else:
                            self.stab_lines[self.stab_tc.create_rectangle(x+1.5*dx,starty_negative,lastx+1.5*dx,endy+starty_negative,
                                                                          fill=self.colour_order[colour_count],
                                                                          outline='',
                                                                          stipple='gray50',
                                                                          width=self.linewidth)]=1
                            label_position[group]=(starty_negative*2+endy)/2.0
                            starty_negative=endy+starty_negative
                    #
                    # Update the colour
                    #
                    colour_count=colour_count+1
                    if colour_count==len(self.colour_order):
                        colour_count=0
            
            lastval=val
            lastpH=pH
        #
        # Put labels on the contributions
        #
        if self.show_grp_contribs.get()==1:
                colour_count=0
                for group in ufgroups:
                    #
                    # Was this group shown?
                    #
                    if group in show_groups:
                        x=self.get_x(pHvalues[-1])
                        y=label_position[group]
                        colour=self.colour_order[colour_count]
                        self.stab_lines[self.stab_tc.create_text(x+1,y,text=group,
                                                           fill=colour,
                                                           anchor='w')]=1
                    #
                    # Update the colour
                    #
                    colour_count=colour_count+1
                    if colour_count==len(self.colour_order):
                        colour_count=0
        #
        #
        #
        lastpH=pHvalues[0]
        lastval=intcurve2[0]
        #
        # Electrostatic contribution
        #
        # THIS IS WRONG!!
        ## count=1
##         for pH in pHvalues[1:]:
##             lastx=self.get_x(lastpH)
##             lasty=pKa_system.get_y_fromstab(lastval,span)
##             val=intcurve2[count]
##             count=count+1
##             x=self.get_x(pH)
##             y=pKa_system.get_y_fromstab(val,span)
##             self.stab_lines[self.stab_tc.create_line(lastx,lasty,float(x),float(y),
##                                                   fill='red',
##                                                   width=self.linewidth)]=1
##             lastval=val
##             lastpH=pH
        #
        # Set focus to the stability window
        #
        self.stab_window.focus_set()
        self.stab_window.lift()
        return

    #
    # ------
    #

    def stability_on_off(self):
        """Open or destroy stability curve window. Called when the stability button is pressed"""
        #
        # Should we open the stability window?
        #
        print 'Opening stab window'
        new=self.stability_var.get()
        if new=='on' and self.old_stab_status!='on':
            self.old_stab_status='on'
            self.stab_window=Toplevel()
            
            #
            #
            # Which contribution should we show
            #
            self.contrib_type=IntVar()
            self.contrib_type.set(1)
            Radiobutton(self.stab_window,
                        text='contributions from pKa shifts',
                        variable=self.contrib_type,value=1,
                        command=self.update_display_calc).grid(row=0,column=1)
            #Radiobutton(self.stab_window,
            #            text='charge-charge contributions',
            #            variable=self.contrib_type,value=2,
            #            command=self.update_display_calc).grid(row=0,column=2)
            
            #
            # Button for showing detailed contributions
            #
            self.show_grp_contribs=IntVar()
            self.show_grp_contribs.set(0)
            grp_contribs=Checkbutton(self.stab_window,text='Show individual residue contributions',
                                     onvalue=1,offvalue=0,variable=self.show_grp_contribs,command=self.update_display_calc)
            grp_contribs.grid(row=0,column=0)
            
            
            #
            # Cutoff for displaying stability
            #
            self.stab_cutoff=DoubleVar()
            self.show_stab_cutoff=Scale(self.stab_window,from_=1,to=50,resolution=1,
                                        orient='horizontal',relief='ridge',
                                        command=self.update_display_calc,variable=self.stab_cutoff,
                                        label='Number of groups to display.',
                                        length=220)
            self.show_stab_cutoff.grid(row=0,column=3)
            self.stab_cutoff.set(10)
            #
            # Do the graphics
            #
            self.stab_window.geometry('+10+430')
            self.stab_window.title('Stability Curve')
            self.stab_tc=Canvas(self.stab_window,bd=5,bg='white',
                                width=self.titwidth,
                                height=self.titheight,
                                scrollregion=(0,0,self.titwidth,self.titheight))
            self.stab_tc.xview("moveto", 0)
            self.stab_tc.yview("moveto", 0)
            self.stab_tc.grid(row=1,column=0,columnspan=4)
            # pH axis
            self.stab_startx,dummy=self.get_xy(self.pHstart,1)
            self.stab_endx,dummy=self.get_xy(self.pHend,1)
            y=pKa_system.get_y_fromstab(0.0,100.0)
            self.stab_starty=160
            self.stab_endy=10
            self.pH_axis(self.stab_tc,self.stab_startx,y,
                         self.stab_endx,y)
            #
            # Button for saving the plot as a postscript file
            #
            def saveplot():
                import tkFileDialog
                filename =tkFileDialog.asksaveasfilename(title='Select filename for postscript file',defaultextension='.ps',
                    initialdir=self.savedir,
                    filetypes=[("Postscript file","*.ps")],parent=self.stab_window)
                if filename:
                    self.stab_tc.postscript(colormode='color',file=filename)
                return
            Button(self.stab_window,text="Save plot",command=saveplot).grid(row=0,column=2)

            #
            # Update
            #
            self.do_stab_curve()
            
        else:
            self.old_stab_status='off'
            self.show_grp_contribs.set(0)
            self.stab_window.destroy()
        return
        
