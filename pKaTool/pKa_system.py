#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
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
from Tkinter import *
import tkFileDialog
import pKa_base, pKa_system_help, pKa_calc, pKa_system_micro, pKa_system_file_menu, CCPS_stab_opt
import pKa_system_data_manipulation
import group_control
import ftir_data

__pKaSystemVersion__=1.2

#
# Geometry helper functions
#

def get_y_fromstab(val,span):
    """Get the y coordinate for plotting the stability curve"""
    zero=10
    graphrange=180
    if span==0:
        span=10
    return (graphrange-val*(graphrange/span))+zero


#
# --------
#


class pKa_system(Frame,pKa_base.pKa_base,pKa_system_help.system_help,
                 pKa_system_help.pKsensitivity,
                 pKa_system_help.decompose,
                 pKa_system_micro.Micro_states,
                 pKa_system_file_menu.file_menu,
                 pKa_system_file_menu.fit_menu,
                 pKa_system_data_manipulation.data_manipulation,
                 CCPS_stab_opt.Optimisation_Analysis):

    def __init__(self,numgroups=None,parent_application=None,data=None,protein=None,field_name=None,update=True):
        #
        # Set up the main window
        #
        # The main window provides choice between the different modes
        # for now: pKa Calculations, pKa System
        #
        self.ID='pKa_system'
        self.font="Times 12 bold"
        self.fg_colour='white'
        self.bg_colour='black'
        self.colour_order=['#646464','#4444ff','red','green','magenta','yellow','orange','grey','magenta']
        self.linewidth=2
        self.names={}
        self.lastdir=None
        self.parent_application=parent_application
        self.init_not_done=1
        self.ID='pKa_system'
        self.protein=protein
        self.field_name=field_name
        #
        # Stability window parameters
        #
        self.stab_window=None
        self.U_control=None
        self.old_stab_status=''
        #
        # Set the pKa calculation parameters
        #
        self.pHstart=0.00001
        self.pHend=20.0
        self.maxcrg=1.0
        self.mincrg=-1.0
        #
        # pKa method
        #
        self.pkamethods={'Boltzmann':pKa_calc.Boltzmann,
                         'Monte Carlo':pKa_calc.Monte_Carlo,
                         'Tanford-Roxby':pKa_calc.Tanford_Roxby,
                         'Monte Carlo (C++)':pKa_calc.Monte_Carlo_CPP,
                         'Boltzmann (C++)':pKa_calc.Boltzmann_CPP}
        #
        # Set data
        #
        self.data={'numgroups':numgroups}
        #
        # All lines drawn
        #
        self.lines={}
        self.state_lines={}
        self.stab_lines={}
        self.to_clear=[]
        self.stab_test_on=None
        #
        # Do the window
        #
        self.do_window()
        if update:
            self.update_pkasystem_curves()
            self.window.update_idletasks()
        #
        # convert data
        #
        if data:
            if data.has_key('groups'):
                self.unpack_all_data(data)
            else:
                self.convert_titration_data(data)
        return

    #
    # ------
    #

    def do_window(self):
        #
        # Main window
        #
        if not self.parent_application:
            Frame.__init__(self)
            self.window=self.master
        else:
            self.window=Toplevel()
        #
        # Title
        #
        #self.window.geometry('+280+380')
        self.window.title("pKa System - Play with titratable groups")
        #
        # Get the size of the screen
        #
        #
        # Text box til top
        #
        label1=Label(self.window, text="Enter number of titratable groups",font=self.font)
        label1.grid(row=0,column=0, sticky=W)
        #
        # Entry field
        #
        self.number_of_groups=IntVar()
        self.numgrps_widget=Entry(self.window,textvariable=self.number_of_groups)
        #
        # If we have a number then insert it
        #
        if self.data['numgroups']:
            self.number_of_groups.set(self.data['numgroups'])
            self.numgrps=self.data['numgroups']
            self.window.destroy()
        else:
            def destroy(event=None):
                self.numgrps=int(self.number_of_groups.get())
                self.window.destroy()
                return
            self.number_of_groups.set(3)
            self.numgrps_widget.grid(row=0,column=1,sticky=W)
            self.numgrps_widget.bind('<Return>',destroy)
            self.window.wait_window(self.window)
        #
        # Done
        #
        self.getgrps()
        return

    #
    # --------------
    #

    def getgrps(self,event=None):
        #
        # Get the number of groups
        #
        import string
        #
        # Open the window for the titration curves
        #
        if not self.parent_application:
            Frame.__init__(self)
            self.window=self.master
            screen_width=self.winfo_screenwidth()
            screen_height=self.winfo_screenheight()
        else:
            self.window=Toplevel()
            screen_width=self.parent_application.winfo_screenwidth()
            screen_height=self.parent_application.winfo_screenheight()

        #
        self.titwin=self.window
        self.titwin.title('Titration curves [native]')
        self.titwin.geometry('+20+%d' %(95+self.numgrps*43))

        #
        # Draw the window with the titration curves
        #
        self.titwidth=1200
        self.titheight=450
        self.tc=Canvas(self.titwin,bd=5,bg='white',width=self.titwidth,
                       height=self.titheight,
                       scrollregion=(0,0,self.titwidth,self.titheight))
        self.tc.xview("moveto", 0)
        self.tc.yview("moveto", 0)
        self.tc.grid(row=0,column=0)
        #
        # Axes
        #
        self.draw_ordinates(self.tc)
        #
        # Open window with the controls
        #
        self.startrow=3
        self.groups={}
        self.win=Toplevel()
        #
        # Create the main pulldown menu
        #
        self.menu=Menu(self.win)
        #
        # File menu
        #
        self.file_menu=Menu(self.menu,tearoff=0)
        self.file_menu.add_command(label='Load system',command=self.load_system)
        self.file_menu.add_command(label='Save system',command=self.save_system)
        if self.parent_application:
            if getattr(self.parent_application,'ID',None):
                if self.parent_application.ID=='EAT':
                    self.file_menu.add_command(label='Save system in EAT & Exit',command=self.send_system_to_EAT)
                elif self.parent_application.ID=='Ekin':
                    self.file_menu.add_command(label='Save system in EAT',command=self.send_system_to_EAT)
        self.file_menu.add_command(label='Load titration curves',command=self.load_curves)
        self.file_menu.add_command(label='Save titration curves',command=self.save_curves)
        self.file_menu.add_command(label='Load titration_DB data',command=self.load_titdb)
        self.file_menu.add_command(label='Load pH activity profile',command=self.load_pH_activity_profile)
        self.file_menu.add_command(label='Load FTIR data',command=self.load_FTIR_data)
        self.file_menu.add_command(label='Print population table',command=self.print_table)
        self.file_menu.add_command(label='Add group',command=self.add_group)
        self.file_menu.add_command(label='Remove exp. titration curve',command=self.remove_exp_curve)
        self.file_menu.add_command(label='Exit',command=self.quit)
        self.menu.add_cascade(label='File',menu=self.file_menu)
        #
        # Command menu
        #
        self.command_menu=Menu(self.menu,tearoff=0)
        self.command_menu.add_command(label='Decompose system',command=self.decompose_system)
        self.command_menu.add_command(label='Sensitivity analysis',command=self.sensitivity_test)
        self.command_menu.add_command(label='Change dielectric constant',command=self.change_dielectric)
        self.command_menu.add_command(label='Activate updating',command=self.activate_callbacks)
        self.command_menu.add_command(label='Deactivate updating',command=self.deactivate_callbacks)
        self.command_menu.add_separator()
        self.command_menu.add_command(label='Copy group to EAT_DB Ekin',command=self.copy_to_Ekin)
        self.menu.add_cascade(label='Command',menu=self.command_menu)
        #
        # Fitting menu
        #
        self.fit_menu=Menu(self.menu,tearoff=0)
        self.fit_menu.add_command(label='Fit system to loaded curves',command=self.fit_system_to_curves)
        self.fit_menu.add_command(label='Fit to pH-activity profile',command=self.fit_to_ph_activity_profile)
        self.fit_menu.add_command(label='Fit to loaded curves and pH-activity profile',command=self.fit_to_curves_and_ph_activity_profile)
        self.fit_menu.add_command(label='Estimate experimental uncertainty',command=self.estimate_experimental_uncertainty)
        self.fit_menu.add_separator()
        self.fit_menu.add_command(label='Fit to FTIR data',command=self.fit_ftir)
        self.fit_menu.add_command(label='Fit to FTIR data and pH-activity profile',command=self.fit_to_ftir_and_ph_activity_profile)
        self.fit_menu.add_separator()
        self.fit_menu.add_command(label='Combinatorial scan',command=self.combinatorial_scan)
        self.fit_menu.add_command(label='Show close parameter sets',command=self.show_close)
        self.fit_menu.add_command(label='Test uniqueness',command=self.test_uniqueness)
        self.fit_menu.add_command(label='Uniqueness scan',command=self.uniqueness_scan)
        self.fit_menu.add_separator()
        self.fit_menu.add_command(label='Evaluate fit',command=self.evaluate_fit)
        self.geom_var=StringVar()
        self.fit_menu.add_command(label="Do geometry optimisation",command=self.start_geom)
        self.fit_menu.add_command(label='Identify number of groups in system',command=self.identify_no_groups)
        self.menu.add_cascade(label='NMR',menu=self.fit_menu)
        #
        # System analysis and optimisation
        #
        self.optana_menu=Menu(self.menu,tearoff=0)
        self.optana_menu.add_command(label='CCPS population/Stability',command=self.stab_and_CCPS_pop)
        self.menu.add_cascade(label='Optimise and Analyse',menu=self.optana_menu)
        #
        # View menu
        #
        self.view_menu=Menu(self.menu,tearoff=0)
        # Show microscopic states
        self.micro_var=IntVar()
        self.micro_var.set(0)
        self.view_menu.add_checkbutton(label='Microscopic titration',
                                        command=self.update_pkasystem_curves,
                                        variable=self.micro_var,onvalue=1,offvalue=0)

        # Show loaded titration curves
        self.display_loaded_curves=IntVar()
        self.display_loaded_curves.set(0)
        self.view_menu.add_checkbutton(label='Loaded titration curves',
                                       command=self.update_pkasystem_curves,
                                       variable=self.display_loaded_curves,onvalue=1,offvalue=0)

        # Show ftir window
        self.show_ftir=IntVar()
        self.show_ftir.set(0)
        self.view_menu.add_checkbutton(label='Show FTIR window',
                                        command=self.update_pkasystem_curves,
                                        variable=self.show_ftir,onvalue=1,offvalue=0)
        #
        # Window for manipulating kcat of microstates
        #
        self.kcat_window_visible=IntVar()
        self.kcat_window_visible.set(0)
        self.view_menu.add_checkbutton(label='kcat of microstates',
        command=self.toggle_kcat_window,variable=self.kcat_window_visible,onvalue=1,offvalue=0)

        self.menu.add_cascade(label='View',menu=self.view_menu)
        #
        # Help menu
        #
        self.help_menu=Menu(self.menu,tearoff=0)
        self.help_menu.add_command(label='About pKaSystem',command=self.about)
        self.menu.add_cascade(label='Help',menu=self.help_menu)
        #
        # Configure the menu
        #
        self.win.config(menu=self.menu)
        #
        self.win.title('Group controls')
        #
        # Place window close to center
        #
        self.win.geometry('+%d+%d' %(screen_width/2-min(600,screen_width/4),screen_height/2-min(500,screen_height/3)))
        #
        # Buttons for each group
        #
        int_ene=1
        for id_num in range(self.numgrps):
            #colour=self.colour_order[id_num%len(self.colour_order)]
            int_ene=1
            self.groups[id_num]=group_control.group_control(self,
                                                            self.startrow+id_num,
                                                            id_num,
                                                            self.numgrps,
                                                            int_ene,
                                                            self.colour_order)
        #
        # Button for controlling stability window
        #
        self.stability_var=StringVar()
        self.stability_var.set('off')
        self.stab_button=Checkbutton(self.win,text='Stability curve: ',
                                     variable=self.stability_var,onvalue='on',
                                     offvalue='off',
                                     command=self.stability_on_off)
        self.stab_button.deselect()
        self.stab_button.grid(row=0,column=0,columnspan=2)
        #
        # Exit button
        #
        self.exit_bt=Button(self.win,text='Quit',command=self.quit)
        self.exit_bt.grid(row=0,column=2,sticky='wens')
        #
        # Snapshot button
        #
        self.snapshot_btn=Button(self.win,text='Snapshot',command=self.snapshot)
        self.snapshot_btn.grid(row=0,column=3,sticky='wens')
        #
        # Window capture button
        #
        self.print_btn=Button(self.win,text='Print2File',command=self.print2file)
        self.print_btn.grid(row=0,column=4,sticky='wens')
        #
        # Clear button
        #
        self.clr_btn=Button(self.win,text='Clear all',command=self.clear)
        self.clr_btn.grid(row=0,column=5,sticky='wens')
        #
        # pHstep slider
        #
        self.pHstep=DoubleVar()
        self.pHstep_sl=Scale(self.win,from_=0.01,to=2.0,resolution=0.01,
                             orient='horizontal',relief='ridge',
                             command=self.update_pkasystem_curves,variable=self.pHstep,
                             label='pHstep')
        self.pHstep_sl.grid(row=0,column=6,sticky='wens')
        self.pHstep.set(0.5)
        #
        # pKa calculation method selector
        #
        self.pkamethod_sel=StringVar()
        self.pkamethod_sel.set('Boltzmann')
        self.pkamethod_button=Menubutton(self.win,textvariable=self.pkamethod_sel,relief=RAISED)
        self.pkamethod_menu=Menu(self.pkamethod_button,tearoff=0)
        self.pkamethod_button['menu']=self.pkamethod_menu
        #
        # Methods
        #
        for method in self.pkamethods.keys():
            self.pkamethod_menu.add_radiobutton(label=method,
                                                variable=self.pkamethod_sel,
                                                value=method,
                                                indicatoron=1,
                                                command=self.update_pkasystem_curves)
        self.pkamethod_button.grid(row=0,column=7,sticky='news')
        #
        # Monte Carlo steps
        #
        self.MCsteps=IntVar()
        self.mcsteps_scale=Scale(self.win,from_=0,to=2500,resolution=100,
                                 orient='horizontal',relief='ridge',
                                 command=self.update_pkasystem_curves,
                                 variable=self.MCsteps,
                                 state=DISABLED,
                                 label='Monte Carlo steps')
        self.MCsteps.set(self.numgrps*100)
        self.mcsteps_scale.grid(row=0,column=8,sticky='news')
        #
        # Button for updating the titration curves
        #
        stab_test=Button(self.win,text='Update curves',command=self.update_pkasystem_curves)
        stab_test.grid(row=0,column=9,sticky='wens')
        #
        # Reposition the window with the titration curves according to the
        # size of the control window
        #
        width,height,xorg,yorg=self.get_geometry(self.win)
        self.titwin.geometry('+%d+%d' %(xorg-5,yorg+height+5))
        #
        # Draw the first curves
        #
        self.titwin.update()
        self.win.update()
        self.init_not_done=None
        self.titwin.update()
        self.win.update()
        self.activate_callbacks()
        #if self.master.update:
        #    self.update_pkasystem_curves()
        #
        # Done
        #
        return
    #
    # -----
    #


    def about(self):
        """Print the About section"""
        import tkMessageBox
        tkMessageBox.showinfo("pKaTool / pKaSystem",
                              'pKaTool version %s\nAuthors: Jens Erik Nielsen & Chresten S¯ndergaard\n\nCopyright (c) Jens Erik Nielsen\nUniversity College Dublin 2003-2007\nAll rigths reserved\nhttp://enzyme.ucd.ie/Science/pKa\n\nPlease remember to cite:\nAnalysing the pH-dependent properties of proteins using pKa calculations\nNielsen JE\nJ Mol Graph Model 2007 Jan;25(5):691-9\n\nIf using the NMR fitting routines please cite:\n\nDetermination of electrostatic interaction energies and protonation state populations in enzyme active sites\nS¯ndergaard CR, McIntosh LP, Pollastri G, Nielsen JE\nJ. Mol. Biol. (in press).' %__pKaSystemVersion__,parent=self.master)
        return


    #
    # --------------------
    #

    def get_geometry(self,widget):
        """Get the geometry of a widget
        Return width,height,xorg,yorg"""
        widget.update_idletasks()
        txt=widget.winfo_geometry()
        width=int(txt.split('x')[0])
        rest=txt.split('x')[1]
        height=int(rest.split('+')[0])
        xorg=int(rest.split('+')[1])
        yorg=int(rest.split('+')[2])
        return width,height,xorg,yorg

    #
    # -------
    #

    def quit_application(self):
        """Quit application"""
        self.win.destroy()
        self.titwin.destroy()
        return

    #
    # --------------------
    #

    def snapshot(self):
        #
        # Preserve the current lines (for making figures)
        #
        x=0
        for line in self.lines.keys():
            x=x+1
            if x==1:
                self.tc.delete(line)
            else:
                self.to_clear.append(line)
            del self.lines[line]
            if x==2:
                x=0
        return

    #
    # --------------------
    #

    def print2file(self):
        #
        # Print Canvas to file
        #
        import sys, os
        if not self.lastdir:
            self.lastdir=os.getcwd()
        filename=tkFileDialog.asksaveasfilename(defaultextension='.ps',
                                                initialdir=self.lastdir,
                                                filetypes=[("Postscript files","*.ps"),("All files","*.*")])
        if filename:
            self.write_psfile(filename)
        else:
            return
        return

    #
    # --------------------
    #

    def write_psfile(self,filename):
        """
        # Dump the Canvas to a postscript file
        """
        import os
        self.lastdir=os.path.split(filename)[0]
        if filename[-3:]!='.ps':
            filename=filename+'.ps'
        self.tc.postscript(colormode='color',file=filename)
        return

    #
    # --------------------
    #

    def clear(self,junk=None):
        #
        # Clear all lines
        #
        for line in self.to_clear:
            self.tc.delete(line)
        return
    #
    # --------------------
    #

    def stability_on_off(self):
        """Open a window for drawing the stability curve"""
        #
        # Should we open the stability window?
        #
        new=self.stability_var.get()
        if new=='on' and self.old_stab_status!='on':
            #
            # Yes, open it
            #
            self.stab_test_on=1
            self.old_stab_status='on'
            self.stab_window=Toplevel()

            #
            self.stabwidth=1000
            self.stabheight=300
            self.stab_window.geometry('%dx%d+10+20' %(self.stabwidth,self.stabheight))
            self.stab_window.title('pH dependence of protein stability')
            self.stab_tc=Canvas(self.stab_window,bd=5,bg='white',width=self.titwidth,height=self.titheight,scrollregion=(0,0,self.titwidth,self.titheight))
            self.stab_tc.xview("moveto", 0)
            self.stab_tc.yview("moveto", 0)
            self.stab_tc.grid(row=1,column=0)
            #
            # Plotting button
            #
            def print_curve(event=None):
                    Button(self.stab_window,command=print_curve).grid(row=0,column=0)
            # pH axis
            self.stab_startx=80
            self.stab_endx=910
            self.stab_starty=160
            self.stab_endy=10
            self.pH_axis(self.stab_tc,self.stab_startx,self.stab_starty,
                         self.stab_endx,self.stab_endy)
            #
            # Controls for unfolded pKa values
            #
            self.U_control=Toplevel()
            self.U_control.title('Controls for Unfolded form')
            self.U_control.geometry('+10+10')
            self.unfolded_groups={}
            for id_num in range(self.numgrps):
                int_ene=1
                #colour=self.colour_order[id_num%len(self.colour_order)]
                self.unfolded_groups[id_num]=group_control.group_control(self,
                                                                         self.startrow+id_num,id_num,
                                                                         self.numgrps,int_ene,self.colour_order,window=self.U_control)
                #
                # If we are displaying real groups then set the intrinsic pKa to the model pKa value
                #
                if self.parent_application:
                    intpKa_folded=self.groups[id_num].modelpK
                else:
                    intpKa_folded=self.groups[id_num].intpka.get()
                #
                self.unfolded_groups[id_num].intpka.set(intpKa_folded)
            #
            #
            #
            row=self.startrow+self.numgrps+1
            self.show_grp_contribs=IntVar()
            self.show_grp_contribs.set(0)
            grp_contribs=Checkbutton(self.U_control,text='Show residue contributions',
                                     onvalue=1,offvalue=0,variable=self.show_grp_contribs,command=self.update_pkasystem_curves)
            grp_contribs.grid(row=row,column=0)
            #
            # Which contribution should we draw
            #
            self.contrib_type=IntVar()
            self.contrib_type.set(1)
            Radiobutton(self.U_control,text='contributions from pKa shifts',variable=self.contrib_type,value=1,command=self.update_pkasystem_curves).grid(row=row,column=1)
            #Radiobutton(self.U_control,text='charge-charge contributions',variable=self.contrib_type,value=2,command=self.update_curves).grid(row=row,column=2)
            #
            # Should we show min and max stabilisation?
            #
            self.show_min_max_stab=IntVar()
            self.show_min_max_stab.set(1)
            Checkbutton(self.U_control,text='Show min and max stabilisation',
                        onvalue=1,
                        offvalue=0,
                        variable=self.show_min_max_stab,
                        command=self.update_pkasystem_curves).grid(row=row,column=3)
            #
            # Update curves
            #
            self.window.update()
            self.U_control.update()
            self.stab_test_on=None
            self.update_pkasystem_curves()
            #
            # Move the windows to sensible positions
            #
            width,height,xorg,yorg=self.get_geometry(self.win)
            self.U_control.geometry('+%d+%d' %(xorg,yorg+height))
            #
            width,height,xorg,yorg=self.get_geometry(self.U_control)
            self.stab_window.geometry('+%d+%d' %(xorg,yorg+height))
            #
            # Activate the callbacks for the unfolded groups
            #
            self.activate_callbacks()
        else:
            self.old_stab_status='off'
            self.stab_window.destroy()
            self.U_control.destroy()
        return
    #
    # --------------------
    #


    def dummy(self,event=None):
        """Dummy callback function"""
        return

    #
    # ----
    #

    def setup_system(self,group_array,X,energies=None):
        """Set up the system of titratable groups from the info in group_array"""
        import string
        #
        # Create the description of the system
        #
        self.names={}
        self.ids={}
        for group in group_array.keys():
            name=':'+string.zfill(group,4)+':'
            if group_array[group].acid_base.get()==1:
                name=name+'ASP'
            else:
                name=name+'ARG'
            self.names[group]=name
            self.ids[name]=group
        #
        # Update experimtal data dictionary to new names...
        #
        if getattr(self,'titration_data',None):
            for old_key in self.titration_data.keys():
                for new_key in self.names:
                    if int(old_key[1:5]) == int(self.names[new_key][1:5]):
                        nk = old_key[0:6]+self.names[new_key][6:]
                        self.titration_data[self.names[new_key]]=self.titration_data[old_key]
                        if not old_key == self.names[new_key]:
                            del self.titration_data[old_key]
        #
        # Add all data
        #
        matrix={}
        X.intene={}
        X.intrinsic_pKa={}
        for group in group_array.keys():
            #
            # Set everything
            #
            name1=self.names[group]
            # Int pKa
            intpka=group_array[group].intpka.get()
            X.intrinsic_pKa[name1]=intpka
            type=group_array[group].acid_base.get()
            #
            # Set the interaction energies
            #
            if not X.intene.has_key(group):
                X.intene[name1]={}
                matrix[name1]={}
                for group2 in group_array[group].intenes.keys():
                    type2=group_array[group2].acid_base.get()
                    name2=self.names[group2]
                    if group_array[group].active.get()==1 and group_array[group2].active.get()==1:
                        if type==type2:
                            X.intene[name1][name2]=group_array[group].intenes[group2].get()
                        else:
                            X.intene[name1][name2]=-group_array[group].intenes[group2].get()
                        if name1!=name2:
                            matrix[name1][name2]=self.E2dist(X.intene[name1][name2],energies)
                    else:
                        X.intene[name1][name2]=0.0
            #
            # We only have part of the interaction energies in each group
            # This is because the interaction energy is stored as a single
            # Tk variable
            #
            for group2 in group_array.keys():
                name2=self.names[group2]
                type2=group_array[group2].acid_base.get()
                if group2!=group:
                    if group_array[group2].intenes.has_key(group):
                        #
                        # Is this group active?
                        #
                        if group_array[group].active.get()==1 and group_array[group2].active.get()==1:
                            if type==type2:
                                X.intene[name1][name2]=group_array[group2].intenes[group].get()
                            else:
                                X.intene[name1][name2]=-group_array[group2].intenes[group].get()
                            #
                            # Matrix of distances
                            #
                            if name1!=name2:
                                matrix[name1][name2]=self.E2dist(X.intene[name1][name2],energies)
                        else:
                            X.intene[name1][name2]=0.0
            #
                else:
                    X.intene[name1][name2]=0.0
                    # Default distance for zero interaction energy
                    if name1!=name2:
                        matrix[name1][name2]=self.E2dist(0.0,energies)
        #
        # All Done
        #
        return matrix

    #
    # -----------------
    #

    def E2dist(self,E,energies=None):
        """
        # convert an electrostatic interaction energy to a distance
        # Units: E(kT), dist: A
        If energies==1, then we do not convert the energy"""
        #
        # Check if we should return energies
        #
        import math
        if energies:
            return E
        #
        # No, return distances
        #
        E=abs(E) # the sign doesn't matter
        if E>0.001:
            #
            # Look in Tynan-Connolly and Nielsen, Protein Science: Re-Designing protein pKa values
            # for details on the formula below
            #
            eps=1.0 # We set eps to 1, and scale distances afterwards
            distance=243.3*math.log(10.0)/(eps*E)
        else:
            distance=1000.0
        return distance

    #
    # ------------------
    #

    def calc_pKas_from_scales(self,group_array):
        """Calculate pKa values for the system"""
        #
        # Fill instance with data
        #
        X=self.pkamethods[self.pkamethod_sel.get()]()
        MCsteps=0
        if self.pkamethod_sel.get()=='Monte Carlo':
            MCsteps=self.MCsteps.get()
            self.mcsteps_scale.configure(state=ACTIVE)
        elif self.pkamethod_sel.get()=='Monte Carlo (C++)':
            MCsteps=200000
        else:
            self.mcsteps_scale.configure(state=DISABLED)
        #
        matrix_dummy=self.setup_system(group_array,X)
        #
        # Set the pKa value variables
        #
        X.groups=X.intrinsic_pKa.keys()
        X.groups.sort()
        #
        # Make a list of experimental pH values to include in calculation
        #
        exp_pHs =[]
        if getattr(self,'titration_data',None):
            for group in self.titration_data.keys():
                for pH in self.titration_data[group].keys():
                    if exp_pHs.count(pH) == 0:
                        exp_pHs.append(pH)
        #
        # also include pH values from loaded ph-activity profile
        #
        if getattr(self,'activity_data',None):
            for pH in self.activity_data.keys():
                if exp_pHs.count(pH) == 0:
                        exp_pHs.append(pH)
        #
        # and also from ftir data
        #
        if getattr(self, 'FTIR_win',None):
            for pH in self.FTIR_win.ftir_data.keys():
                if exp_pHs.count(pH) ==0:
                    exp_pHs.append(pH)
        #
        # Include the effect of non-system groups?
        #
        if hasattr(self,'non_system_groups'):
            if self.non_system_groups:
                X.non_system_groups={}
                for group_id in self.non_system_groups.keys():
                    name=self.names[group_id]
                    X.non_system_groups[name]=self.non_system_groups[group_id].copy()
        #
        # Get the pKa values
        #
        pKa_values,prot_states=X._calc_pKas(mcsteps=MCsteps,
                                            phstep=self.pHstep.get(),
                                            phstart=self.pHstart,
                                            phend=self.pHend,
                                            exp_pHs=exp_pHs,
                                            verbose=1)
        return X,pKa_values,prot_states

    #
    # -----------------
    #

    def update_scales(self,junk=None,draw=1,doit=None):
        """Update the scale widgets when the user moves a dial"""
        #
        # Folded (normal) groups
        #
        for group in self.groups.keys():
            self.groups[group].update_scales()
        #
        # Update the unfolded scales if theyr're active
        #
        if self.stability_var.get()=='on':
            for group in self.unfolded_groups.keys():
                self.unfolded_groups[group].update_scales()
        #
        # Redraw everything
        #
        self.update_pkasystem_curves(junk,draw,doit)
        return

    #
    # -----
    #

    def update_scales_from_fit(self,junk=None,draw=1,doit=None):
        #
        # update group scales from fitter
        #
        for group in self.groups:
            self.groups[group].update_scales_from_fit()
        self.update_pkasystem_curves(junk,draw,doit)
        return

    #
    # -----
    #

    def update_pkasystem_curves(self,junk=None,draw=1,doit=None):
        """Update all curves"""
        if self.init_not_done:
            return
        if self.stab_test_on and doit==None:
            return
        #
        # Redraw the curves
        #
        import string, pKarun
        PKana=pKarun.pKa_general.pKanalyse()
        #
        # Calculate pKa values for the folded form
        #

        X,pKa_values,prot_states=self.calc_pKas_from_scales(self.groups)
        self.pKa_calc_instance=X
        if not draw:
            return X
        #
        # Set the pKa values
        #
        for group in pKa_values.keys():
            self.groups[self.ids[group]].update_group_control()
            self.groups[self.ids[group]].pkavalue.set("%4.1f" %pKa_values[group])
            #
            # Set the HH fit
            #
            solution,sq=PKana.fit_to_henderson(X.prot_states[group])
            try:
                self.groups[self.ids[group]].HHfit.set('%5.2f (%4.2f / %3.2f)' %(abs(float(solution[1])),abs(float(solution[0])),float(sq)))
            except:
                self.groups[self.ids[group]].HHfit.set('HH-fit error')
        #
        # Delete all lines from last round
        #
        for line in self.lines.keys():
            self.tc.delete(line)
            del self.lines[line]

        # Draw the titration curves
        self.titration_curves={}
        groups=pKa_values.keys()
        groups.sort()
        group_count=0
        colour_map = {}
        for group in groups:
            #
            # Store everything in self.titration_curves
            #
            self.titration_curves[group]=X.prot_states[group].copy()
            #
            # Is this group active?
            #
            if self.groups[group_count].active.get()==0:
                group_count=group_count+1
                continue
            #
            # Yes
            #
            style=self.groups[group_count].style.get()
            lastpH=X.pHvalues[0]
            lastcrg=X.prot_states[group][lastpH]
            colour=self.colour_order[group_count%len(self.colour_order)]
            colour_map[group] = colour

            for pH in X.pHvalues[1:]:
                lastx,lasty=self.get_xy(lastpH,lastcrg)
                crg=X.prot_states[group][pH]
                x,y=self.get_xy(pH,crg)

                if style==1:
                    self.lines[(self.tc.create_line(lastx,lasty,float(x),float(y),
                                                    fill=colour,
                                                    width=self.linewidth))]=1
                else:
                    self.lines[(self.tc.create_line(lastx,lasty,float(x),float(y),
                                                    fill=colour,
                                                    width=self.linewidth,
                                                    dash=(1,2)))]=1
                lastcrg=crg
                lastpH=pH



            #
            # Update the counter for colours
            #
            group_count=group_count+1
        #
        # Should we draw the microscopic states?
        #
        if self.micro_var.get()==1:
            self.update_microstates(X)
        else:
            self.close_state_win()

        #
        # Should we draw the stabilty curves?
        #
        stab_status=self.stability_var.get()
        if stab_status=='on':
            self.stability=self.do_stab_curve(X)
        #
        # Should we display loaded titration curves?
        #
        try:
            print 'titration_data', self.titration_data
        except:
            print 'no titration_data'
        if self.display_loaded_curves.get()==1:
            if not getattr(self,'titration_data',None):
                import tkMessageBox
                tkMessageBox.showwarning('No titration curves loaded',
                                         'Load titration curves first')
                self.display_loaded_curves.set(0)
            else:
                for group in self.titration_data.keys():
                    phvals=self.titration_data[group].keys()
                    phvals.sort()
                    for ph in phvals:
                        crg=self.titration_data[group][ph]
                        x,y=self.get_xy(ph,crg)
                        try:
                            f = colour_map[group]
                        except:
                            f = 'yellow'

                        handle=self.tc.create_oval(x-2,y-2,x+2,y+2,fill=f)

                        self.lines[handle]=1
        #
        # Is there an FTIR model to update?
        #
        if self.show_ftir.get() == 1:
            if getattr(self, 'FTIR_win',None):
                self.FTIR_win.draw_fit()
            else:
                self.FTIR_win = ftir_data.FTIR_data(self)
                self.FTIR_win.draw_fit()
        #
        # Other callbacks?
        #
        self.check_other_callbacks()
        return X

    #
    # -------
    #

    def check_other_callbacks(self):
        """self.callbacks holds a list of funcions that should be called"""
        if not hasattr(self,'callbacks'):
            self.callbacks=[]
        for callback in self.callbacks:
            callback()
        return

    def add_callback(self,function):
        """Add a callback function"""
        self.check_other_callbacks()
        add=1
        for callback in self.callbacks:
            if function==callback:
                add=None
                break
        if add:
            self.callbacks.append(function)
        self.check_other_callbacks()
        return

    #
    # -------------
    #

    def do_stab_curve(self,X):
        """ Calculate the stability curve"""
        #
        # Make sure that the acid/base info for the unfolded form is the same
        # as for the folded form
        #
        for group in self.unfolded_groups.keys():
            acid_base=self.groups[group].acid_base.get()
            self.unfolded_groups[group].acid_base.set(acid_base)
        #
        # Calculate pKa values for the unfolded form
        #
        UF,ufpKa_values,UF_prot_states=self.calc_pKas_from_scales(self.unfolded_groups)
        for group in ufpKa_values.keys():
            self.unfolded_groups[self.ids[group]].pkavalue.set("%4.1f" %ufpKa_values[group])
        #
        # Get all the interaction energies
        #
        ufmatrix=self.setup_system(group_array=self.unfolded_groups,X=UF,energies=1)
        matrix=self.setup_system(group_array=self.groups,X=X,energies=1)
        #
        # Integrate
        #
        integral=0.0
        intcurve=[]
        intcurve2=[]
        dpH=abs(X.pHvalues[0]-X.pHvalues[1])
        min_val=99999
        max_val=-9999
        #
        # Specify constants
        #
        k=1.3806503E-23
        T=298.15
        Na=6.02214199E23
        #factor=k*T*Na/1000.0
        # No don't do it
        factor=1
        import math
        ln10=math.log(10)
        #
        # Loop over all pH values
        #
        stability={}
        for pH in X.pHvalues:
            intcurve.append(integral)
            stability[pH]=integral #Dictionary that will be passed back
            #
            # Update min and max
            #
            if integral<min_val:
                min_val=integral
            if integral>max_val:
                max_val=integral
            #
            # Calculate total stability
            #
            for group in ufpKa_values.keys():
                integral=integral+ln10*dpH*(X.prot_states[group][pH]-UF.prot_states[group][pH])*factor
            #
            # Calculate the electrostatic interaction
            #
            integral2=0
            for group in matrix.keys():
                for group2 in matrix.keys():
                    #
                    # Get the interaction between this group and the other group
                    #
                    g1_id=self.ids[group]
                    g2_id=self.ids[group2]
                    if self.groups[g1_id].active.get()==1 and self.groups[g2_id].active.get()==1 and group!=group2:
                        integral2=integral2+abs(X.prot_states[group][pH])*abs(X.prot_states[group2][pH])*matrix[group][group2]/2.0*factor
                        # Subtract the interaction in the unfolded state
                        integral2=integral2-abs(UF.prot_states[group][pH])*abs(UF.prot_states[group2][pH])*ufmatrix[group][group2]/2.0*factor
            #
            # Update min and max
            #
            if integral2<min_val:
                min_val=integral2
            if integral2>max_val:
                max_val=integral2
            intcurve2.append(integral2)
        max_stabilisation=max_val
        min_stabilisation=min_val
        #
        # Plot the whole thing
        #
        lastpH=X.pHvalues[0]
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
        x_axis=self.get_x(X.pHvalues[0])-20
        y_axis=get_y_fromstab(min_val,span)
        endy=get_y_fromstab(max_val,span)
        self.stab_lines[canvas.create_line(x_axis,max([160,y_axis]),
                                                  x_axis,endy-10,fill='black',
                                                  width=self.linewidth)]=1
        self.stab_lines[canvas.create_text(x_axis+10,endy-35,text='delta G of folding (kT)',fill='black',anchor='w')]=1
        #self.stab_lines[canvas.create_text(x_axis+10,endy-55,text='Net Electrostatic interactions (folded - unfolded) (kT)',
        #                                   fill='red',
        #                                   anchor='w')]=1
        #
        # Tick marks and tick labels
        #
        for tickval in range(int(min_val*100),int(max_val*100),int(max([(span*100.0)/5.0,1.0]))):
            y=get_y_fromstab(tickval/100.0,span)
            self.stab_lines[canvas.create_line(x_axis,
                                                y,x_axis-5,y,
                                               fill='black',width=self.linewidth)]=1
            self.stab_lines[canvas.create_text(x_axis-25,y,text='%5.2f' %(
                float(tickval)/100.0),fill='black')]=1
        #
        # Draw the stability lines
        #
        count=1
        summed_contributions={}
        label_position={}
        for pH in X.pHvalues[1:]:
            lastx=self.get_x(lastpH)
            lasty=get_y_fromstab(lastval,span)
            val=intcurve[count]

            x=self.get_x(pH)
            y=get_y_fromstab(val,span)
            self.stab_lines[self.stab_tc.create_line(lastx,lasty,float(x),float(y),
                                                  fill='black',
                                                  width=self.linewidth)]=1
            #
            # Outline the contribution of each group
            #
            if self.show_grp_contribs.get()==1:
                colour_count=0
                null_y=get_y_fromstab(0.0,span)
                starty_positive=null_y
                starty_negative=null_y
                ufgroups=ufpKa_values.keys()
                ufgroups.sort()
                for group in ufgroups:
                    #
                    # Is this group active?
                    #
                    g1_id=self.ids[group]
                    if self.groups[g1_id].active.get()==1:
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
                            endy=get_y_fromstab(dpH*ln10*(X.prot_states[group][pH]-UF.prot_states[group][pH])*factor,span)-null_y
                            summed_contributions[group]=summed_contributions[group]+endy
                        else:
                            #
                            # Otherwise the stability contribution from charge-charge interactions
                            #
                            stab=0.0
                            for group2 in matrix.keys():
                                #
                                # Get the interaction between this group and the other group
                                #
                                g2_id=self.ids[group2]
                                if self.groups[g1_id].active.get()==1 and self.groups[g2_id].active.get()==1 and group!=group2:
                                    stab=stab+abs(X.prot_states[group][pH])*abs(X.prot_states[group2][pH])*matrix[group][group2]/2.0*factor
                                    # Subtract the interaction in the unfolded state
                                    stab=stab-abs(UF.prot_states[group][pH])*abs(UF.prot_states[group2][pH])*ufmatrix[group][group2]/2.0*factor
                            endy=get_y_fromstab(stab,span)-null_y
                            summed_contributions[group]=endy
                        #
                        # Add it
                        #

                        #
                        # Draw the box
                        #
                        endy=summed_contributions[group]
                        if endy>0:
                            self.stab_lines[self.stab_tc.create_rectangle(x+1.5*dx,starty_positive,lastx+1.5*dx,endy+starty_positive,
                                                                          fill=self.colour_order[colour_count],
                                                                          outline=self.colour_order[colour_count],
                                                                          stipple='gray50',
                                                                          width=self.linewidth)]=1
                            label_position[group]=(starty_positive*2+endy)/2.0
                            starty_positive=endy+starty_positive

                        else:
                            self.stab_lines[self.stab_tc.create_rectangle(x+1.5*dx,starty_negative,lastx+1.5*dx,endy+starty_negative,
                                                                        fill=self.colour_order[colour_count],
                                                                          outline=self.colour_order[colour_count],
                                                                          stipple='gray50',
                                                                          width=self.linewidth)]=1
                            label_position[group]=(starty_negative*2+endy)/2.0
                            starty_negative=endy+starty_negative

                    colour_count=colour_count+1
                    if colour_count==len(self.colour_order):
                        colour_count=0
            #
            # Continue
            #
            lastval=val
            lastpH=pH
            count=count+1
        #
        # Put labels on the contributions
        #
        if self.show_grp_contribs.get()==1:
                colour_count=0
                for group in ufgroups:
                    #
                    # Is this group active?
                    #
                    g1_id=self.ids[group]
                    if self.groups[g1_id].active.get()==1:
                        x=self.get_x(X.pHvalues[-1])
                        y=label_position[group]
                        colour=self.colour_order[colour_count]
                        self.stab_lines[canvas.create_text(x+50,y,text=group,
                                                           fill=colour,
                                                           anchor='w')]=1
                    #
                    # Update colours
                    #
                    colour_count=colour_count+1
                    if colour_count==len(self.colour_order):
                        colour_count=0
        #
        # Now for the electrostatic contribution
        #
        ## count=1
##         lastpH=X.pHvalues[0]
##         lastval=intcurve2[0]
##         count=1
##         for pH in X.pHvalues[1:]:
##             lastx=self.get_x(lastpH)
##             lasty=get_y_fromstab(lastval,span)
##             val=intcurve2[count]

##             x=self.get_x(pH)
##             y=get_y_fromstab(val,span)
##             self.stab_lines[self.stab_tc.create_line(lastx,lasty,float(x),float(y),
##                                                   fill='red',
##                                                   width=self.linewidth)]=1
##             lastval=val
##             lastpH=pH
##             count=count+1
        #
        # Put in labels for min and max stabilisation
        #
        if self.show_min_max_stab.get()==1:
            obj1=canvas.create_text(850,150,text='MAX destab: %5.2f kT' %max_stabilisation,fill='red',anchor='w')
            obj2=canvas.create_text(850,180,text='MAX stab: %5.2f kT' %min_stabilisation,fill='blue',anchor='w')
            self.stab_lines[obj1]=1
            self.stab_lines[obj2]=1
        return stability

    #
    # ---------------
    #

    def start_geom(self):
        """Start geom opt"""
        import pKa_calc
        X=pKa_calc.Boltzmann()
        distance_matrix=self.setup_system(self.groups,X)
        import dist_geom
        GM=dist_geom.distance_optimisation(distance_matrix,self.titration_curves)
        return

    #
    # ----
    #

    def do_geom(self):
         #
         # Do geometry optimisation
         #
         # Update distances
         #
         import pKa_calc
         X=pKa_calc.Boltzmann()
         distance_matrix=self.setup_system(self.groups,X)
         self.MD.set_eqdists(distance_matrix)

         #diff=self.MD.EM(1)
         #
         # Delete old ovals
         #
         for oval in self.oval.keys():
             self.geom_tc.delete(oval)
             del self.oval[oval]
         #
         # Plot positions
         #
         group_count=0
         groups=self.MD.atoms.keys()
         groups.sort()
         for grp in groups:
             pos=self.MD.atoms[grp]['pos']
             x=pos[0]
             y=pos[1]
             z=pos[2]
             self.oval[self.geom_tc.create_oval(x-5,y-5,x+5,y+5,fill=self.colour_order[group_count])]=1
             group_count=group_count+1
             self.oval[self.geom_tc.create_text(10,10,anchor='nw',text='Sum of unsatisfied dists: %5.3f' %(diff))]=1
         self.geom_window.after(100,self.start_geom)
         return

    #
    # ----------------
    #

    def copy_to_Ekin(self):
        """Copy a titration curve or a population curve to the Ekin facility of EAT_DB"""
        try:
            import os,sys
            import EAT_DB.Ekin
        except:
            import tkMessageBox
            tkMessageBox.showwarning('Cannot find PEAT',
                                     'Cannot find PEAT_DB\nMake sure that you download PEAT from\nhttp://enzyme.ucd.ie/PEAT')
            return
        #
        # Pick a group
        #
        self.pick_group=Toplevel()
        self.pick_group.title('Pick a group')
        self.pick_group.geometry('+200+200')
        self.group_picked=IntVar()
        count=0
        groups=self.groups.keys()
        groups.sort()
        for group in groups:
            Radiobutton(self.pick_group,text='%d:%s' %(group,self.groups[group].name.get()),
                        variable=self.group_picked,
                        value=count).grid(row=count,column=0)
            count=count+1
        self.group_picked.set(groups[0])
        Button(self.pick_group,text='Copy group',command=self.copy_group).grid(row=count,column=0)
        Button(self.pick_group,text='Cancel',command=self.cancel_copy_group).grid(row=count,column=1)
        return

    #
    # ----
    #

    def copy_group(self,event=None):
        """Get the titration curve and send it to Ekin"""
        group=None
        for id in self.ids.keys():
            if self.ids[id]==self.group_picked.get():
                group=id
                break
        if not group:
            raise 'Something very odd happended in copy_group'
        #
        # Get the data and reformat it
        #
        data=self.titration_curves[group].copy()
        del data['pKa']
        new_data={}
        new_data[0]={}
        new_data[1]={}
        count=0
        pHs=self.titration_curves[group].keys()
        pHs.sort()
        for pH in pHs:
            new_data[0][count]=pH
            new_data[1][count]=self.titration_curves[group][pH]
            count=count+1
        #
        # Open Ekin, and load the data
        #
        import os,sys
        import EAT_DB.Ekin
        EK=EAT_DB.Ekin.Ekin(parent=self)
        EK.pass_data(new_data)
        #
        # Destroy the little window
        #
        self.pick_group.destroy()
        return

    #
    # ----
    #

    def cancel_copy_group(self,event=None):
        """Cancel copy group to Ekin"""
        self.pick_group.destroy()
        return



#
# -----------------
#

if __name__=='__main__':
    import sys
    if len(sys.argv)==2:
        numgroups=int(sys.argv[1])
        pKa_system(numgroups).mainloop()
    else:
        pKa_system().mainloop()
