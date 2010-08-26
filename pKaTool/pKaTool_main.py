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
import tkMessageBox, tkFileDialog


import os, sys, time

import pKa_base, pKaIO, pKa_system, pKacalc_interface, pI, pKaTool_stability, pKaTool_utility
import pKarun.pKarun_main as pKarun
from pKarun.pKa_utility_functions import *


#
# Version identifier
#
__pKaToolVersion__=1.05

class Main(Frame,pKaTool_stability.pKaTool_stability,pKaTool_utility.pKaTool_utility,pKa_base.pKa_base,pKacalc_interface.Main,pI.pI):

    def __init__(self):
        #
        # Set up the main window
        #
        # The main window provides choice between the different modes
        # for now: pKa Calculations, pKa System
        #
        self.font="Times 12 bold"
        self.fg_colour='white'
        self.bg_colour='black'
        #
        # Open the window
        #
        Frame.__init__(self)
        self.master.title("pKaTool")
        #
        # Place window close to center
        #
        screen_width=self.winfo_screenwidth()
        screen_height=self.winfo_screenheight()
        self.master.geometry('+%d+%d' %(screen_width/2-min(600,screen_width/4),screen_height/2-min(500,screen_height/3)))
        #
        # Stability curve parameters
        #
        self.stab_window=None
        self.U_control=None
        self.old_stab_status=''
        self.stab_lines={}
        #
        # Did we get a pdb file with the commandline?
        #
        self.pdbfile=None
        #
        # A few more variables
        #
        import os
        self.savedir=os.getcwd()
        self.num_residues=0
        self.pkavals={}
        self.pHstart=0.0
        self.pHend=20.0
        self.pHstep=0.1
        self.label_count=0
        self.labels=[]
        for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            self.labels.append(letter)
        for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            for letter2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                self.labels.append(letter+letter2)
        self.pdbfilename=None
        self.pdbfilename_var=StringVar()
        self.pdbfilename_var.set('PDB file: %s' %str(self.pdbfile))
        self.pka_status='N/A'
        self.pka_status_var=StringVar()
        self.pka_status_var.set('pKa calc status: %s' %self.pka_status)
        #
        self.current_selected_calc=None
        self.groups=[]
        self.calcs={}
        #
        # Calculation engine
        #
        self.calc_engine=StringVar()
        self.calc_engine.set('pKaTool')
        self.pdb2pka_path=StringVar()
        self.pdb2pka_path.set("")
        #
        # Go setup the window
        #
        self.analysis_mode()
        return
    #
    # -----
    #

    def quit(self):
        self.master.destroy()
        os._exit(0)
        return

    #
    # ------
    #

    def set_position(self,window):
        """Set the position of a new window to be on top of self.master"""
        geom=self.master.winfo_geometry()
        geom=geom.split('+')
        geom='+%s+%s' %(geom[-2],geom[-1])
        window.geometry(geom)
        return

    #
    # ------
    #

    def pulldown(self):
        """Create the pulldown menu"""
        #
        # Create the main pulldown menu
        #
        self.menu=Menu(self.master)
        #
        # File menu
        #
        self.file_menu=Menu(self.menu,tearoff=0)
        self.file_menu.add_command(label='Load PDB file',command=self.load_pdb)
        self.file_menu.add_command(label='Load calculation(s)',command=self.load_titration)
        self.file_menu.add_command(label='Load calculation(s) from pKD server',command=self.load_frompKD)
        self.file_menu.add_command(label='Save calculations(s)',command=self.save_titration)
        self.file_menu.add_command(label='Subtract calculation',command=self.subtract_titration)
        self.file_menu.add_command(label='Quit',command=self.quit)
        self.menu.add_cascade(label='File',menu=self.file_menu)
        #
        # Calculation menu
        #
        self.calculation_menu=Menu(self.menu,tearoff=0)
        self.calculation_menu.add_command(label='Initialise',command=self.init_calc)
        self.calculation_menu.add_command(label='Desolvation energies',state=DISABLED)
        self.calculation_menu.add_command(label='Background interaction energies',state=DISABLED)
        self.calculation_menu.add_command(label='Interaction energy  matrix',command=self.calculate_matrix_dialog)
        self.calculation_menu.add_command(label='Perform full pKa calculation',command=self.do_pka_calculation)

        self.calculation_menu.add_command(label='Calculate titration curves',state=DISABLED)
        self.calculation_menu.add_separator()
        self.calculation_menu.add_command(label='Set titration curve calculation method',command=self.set_titcurv_method,state=DISABLED)
        self.calculation_menu.add_command(label='Specify paths',state=DISABLED)
        self.menu.add_cascade(label='pKa calculation',menu=self.calculation_menu)
        
        #
        # Select menu
        #
        self.select_menu=Menu(self.menu,tearoff=0)
        self.menu.add_cascade(label='Select',menu=self.select_menu)
        #
        self.display_menu=Menu(self.menu,tearoff=0)
        self.menu.add_cascade(label='Display titration curves',menu=self.display_menu)
        #
        # Analyse menu
        #
        self.analyse_menu=Menu(self.menu,tearoff=0)
        self.analyse_menu.add_command(label='Show table with group info',command=self.show_details_table)
        self.analyse_menu.add_command(label='Multiple calculation analysis',command=self.mult_calc_anal)
        self.analyse_menu.add_command(label='Calculate pI',command=self.calc_pI)
        self.menu.add_cascade(label='Analyse',menu=self.analyse_menu)
        #
        # Help menu
        #
        self.help_menu=Menu(self.menu,tearoff=0)
        self.help_menu.add_command(label='About pKaTool',command=self.about)
        self.menu.add_cascade(label='Help',menu=self.help_menu)
        #
        self.master.config(menu=self.menu)
        #
        # Variable for holding the selected calculation
        #
        self.calc_disp={}
        self.calc_select=StringVar()
        #
        # The currently selected group
        #
        self.current_selected_calc=None
 
        return

    #
    # -----
    #

    def mult_calc_anal(self,event=None):
        """Open the multiple calc analysis window"""
        import multiple_calc_analysis
        multiple_calc_analysis.multiple_calc_analysis(self)
        return

    #
    # -----
    #

    def about(self):
        """Print the About section"""
        tkMessageBox.showinfo("pKaTool",
                              'pKaTool version %s\nAuthors: Jens Erik Nielsen & Chresten Soendergaard\nCopyright (c) Jens Erik Nielsen\nUniversity College Dublin 2003-2007\nAll rigths reserved\nhttp://enzyme.ucd.ie/Science/pKa\n\nPlease remember to cite:\n\nAnalysing the pH-dependent properties of proteins using pKa calculations\nNielsen JE\nJ Mol Graph Model.2007 Jan;25(5):691-9' %__pKaToolVersion__,parent=self.master)
        return

    #
    # -----
    #

    def add_to_pulldown(self,pdbfilename,label):
        """Add a new titration to the menus"""
        self.calc_disp[pdbfilename]=IntVar()
        self.display_menu.add_checkbutton(label=label+':'+pdbfilename,
                                          variable=self.calc_disp[pdbfilename],
                                          onvalue=1,offvalue=0,
                                          command=self.update_display_calc)
        #
        # Now the select menu
        #
        self.select_menu.add_radiobutton(label=label+':'+pdbfilename,
                                         variable=self.calc_select,
                                         value=pdbfilename,
                                         indicatoron=1,
                                         command=self.update_selected_calc)
        return

        
    
    #
    # -----
    #

    def analysis_mode(self):
        """Start the analysis mode of pKaTool. this involves setting up all the gadgets etc.
        The new form of the analysis mode does not need a completed pKa calculation - or a PDB file"""
        #
        # Add to the display and select menus, and select this calculation for both
        #
        self.pulldown()

        #
        # Expand window 
        #
        row=0
        Label(self.master,textvariable=self.pdbfilename_var).grid(row=row,column=0)
        Label(self.master,textvariable=self.pka_status_var).grid(row=row,column=3)
        row=row+1
        Label(self.master,text='Groups in file (%d titratable groups/%d amino acids)' %(len(self.pkavals.keys()),
                                                                                        self.num_residues)).grid(row=row,column=0,columnspan=4,sticky='nws')
        #
        # Put a listbox in for selecting groups
        #
        Label(self.master,text='Selected groups').grid(row=row,column=4)
        Label(self.master,text='Group detail').grid(row=row,column=6)
        #
        # Add the boxes themselves and their scrollbars
        #
        rowspan=8
        row=row+1
        yscrollbar=Scrollbar(self.master,orient='vertical',width=10)
        yscrollbar.grid(row=row,column=2,rowspan=rowspan,sticky='nws')
        #
        height=10
        self.group_disp=Listbox(self.master,bg='white',
                            fg='black',
                            height=height,width=15,
                            yscrollcommand= yscrollbar.set,
                            selectmode=EXTENDED)
        self.group_disp.grid(row=row,column=0,columnspan=2,rowspan=rowspan,sticky='news')
        yscrollbar.config(command=self.group_disp.yview)
        #
        # Buttons for adding and remove groups to the selection
        #
        Button(self.master,text='Add -->',command=self.select_group).grid(row=row,column=3)
        Button(self.master,text='<-- Remove',command=self.remove_group).grid(row=row+1,column=3)
        Button(self.master,text='Show selected subsystem',command=self.show_system).grid(row=row+2,column=3)
        #
        Button(self.master,text='pKaSystem',command=self.goto_pKa_system).grid(row=row+3,column=3)
        Button(self.master,text='Quit',command=self.quit).grid(row=row+4,column=3)
                
        #
        # List containing selected groups
        #
        yscrollbar2=Scrollbar(self.master,orient='vertical',width=10)
        yscrollbar2.grid(row=row,column=5,rowspan=rowspan,sticky='nws')
        self.selected=Listbox(self.master,bg='white',
                            fg='black',
                            height=height,width=30,
                            yscrollcommand= yscrollbar2.set,
                            selectmode=EXTENDED)
        self.selected.grid(row=row,column=4,columnspan=1,rowspan=rowspan,sticky='NWS')
        yscrollbar2.config(command=self.selected.yview)

        #
        # Text area for details on groups
        #
        column=5
        column=column+1
        self.details=Text(width=30,height=15,bg='white')
        self.details.grid(row=row,column=column,rowspan=rowspan)
        row=row+rowspan
        # -------------------------
        #
        # Fields for selecting groups
        #
        row=row+1
        column=0
        Label(self.master,text='Get pKa of').grid(row=row,column=column)
        column=column+1
        choices=['Acids','Bases','Asp + Glu','Lys + Arg','HIS','ARG','CYS','ASP','GLU','TYR','LYS','Any group']
        self.type=StringVar()
        self.type.set('Any group')
        # the button
        self.t_button=Menubutton(self.master,textvariable=self.type,relief=RAISED)
        # the menu
        self.t_menu=Menu(self.t_button,tearoff=0)
        self.t_button['menu']=self.t_menu
        for type in choices:
            self.t_menu.add_radiobutton(label=type,variable=self.type,value=type,indicatoron=1)
        self.t_button.grid(row=row,column=column,sticky='news')
        column=column+2
        #
        # qualifiers
        #
        self.qualifier=StringVar()
        list=['raised by more than','lowered by more than','greater than','smaller than','equal to']
        button=self.menubutton_list(window=self.master,
                                    variable=self.qualifier,
                                    list=list,
                                    default='raised by more than',
                                    indicatoron=0)
        button.grid(row=row,column=column,sticky='news')
        column=column+1
        #
        # Entry value
        #
        self.amount=DoubleVar()
        self.amount.set(2.0)
        Entry(self.master,textvariable=self.amount).grid(row=row,column=column,sticky='news')
        column=column+1
        Button(self.master,text='Search',command=self.find_groups).grid(row=row,column=column,sticky='news')
        # -----------------------------------------------------------------
        #
        # Add button for selecting odd titration curves
        #
        row=row+1
        column=0
        Label(self.master,text='Get non-Henderson-Hasselbalch groups').grid(row=row,column=column,sticky='news',columnspan=2)
        column=column+3
        Label(self.master,text='Cutoff (slope of HH fit)').grid(row=row,column=column)
        column=column+1
        self.irr_cutoff=DoubleVar()
        self.irr_cutoff.set(0.8)
        Entry(self.master,textvariable=self.irr_cutoff).grid(row=row,column=column,sticky='news')
        column=column+1
        Button(self.master,text='Search',command=self.select_irregular_titcurvs).grid(row=row,column=column,sticky='news')
        #
        # ---------------------------------------------------------------------------------------------
        # Button for finding strongly coupled groups
        #
        column=0
        row=row+1
        Label(self.master,text='Add strongly interacting groups').grid(row=row,column=column,columnspan=2,sticky='news')
        #
        column=column+3
        Label(self.master,text='Cutoff (kT/e)').grid(row=row,column=column)
        column=column+1
        self.int_cutoff=DoubleVar()
        self.int_cutoff.set(1.0)
        Entry(self.master,textvariable=self.int_cutoff).grid(row=row,column=column,sticky='news')
        column=column+1
        Button(self.master,text='Add groups',command=self.select_strongly_interacting_groups).grid(row=row,column=column,sticky='news')
        #
        # Checkbox for the stability window
        #
        column=column+1
        self.stability_var=StringVar()
        self.stability_var.set('off')
        self.stab_button=Checkbutton(self.master,text='Show stability profile',
                                     variable=self.stability_var,onvalue='on',
                                     offvalue='off',
                                     command=self.stability_on_off)
        self.stab_button.deselect()
        self.stab_button.grid(row=row-2,column=column)
        #
        # Set parameters for drawing titration curves
        #
        self.mincrg=-1
        self.maxcrg=1
        self.linewidth=2
        #
        # Bottom window for showing titration curves
        #
        row=row+1
        self.titwidth=1200
        self.titheight=400
        self.tc=Canvas(self.master,bd=5,bg='white',width=self.titwidth,
                       height=self.titheight,scrollregion=(0,0,self.titwidth,self.titheight))
        self.tc.xview("moveto", 0)
        self.tc.yview("moveto", 0)
        self.tc.grid(row=row,column=0,columnspan=7)
        #
        # Draw the axes
        #
        self.draw_ordinates(self.tc)
        #
        # Update all groups
        #
        self.update_all_groups()
        #
        # Clear the selection
        #
        self.selected_groups=[]
        #
        # Bind left-click to display of details
        #
        self.group_disp.bind('<ButtonRelease-1>',self.display_details)
        self.selected.bind('<ButtonRelease-1>',self.display_details_sel)
        #
        # Ctrl+A for selecting All
        #
        self.master.bind('<Control-a>',self.select_all)
        return

    #
    # --------
    #

    def update_all_groups(self,event=None):
        """Update the listbox showing the groups for the selected calculation"""
        #
        # Add the groups to the listbox
        #
        self.group_disp.delete(0, END)
        for group in self.groups:
            self.group_disp.insert(END,group)
        return
        
    #
    # --------
    #

    def update_selected_calc(self,event=None):
        """Update the set of pKa values we're working with"""
        #
        # First store the selected groups from the old selection
        #
        if self.current_selected_calc:
            self.calcs[self.current_selected_calc]['selected_groups']=self.selected_groups[:]
        #
        # Then load the new stuff
        #
        new_calc=self.calc_select.get()
        self.X=self.calcs[new_calc]['instance']
        self.titration_curves=self.calcs[new_calc]['titcurv']
        self.pkavals=self.calcs[new_calc]['pkavals'].copy()
        if self.X:
            self.matrix=self.X.read_matrix()
        #
        self.groups=self.titration_curves.keys()
        self.groups.sort()
        self.current_selected_calc=new_calc
        #
        # Update the selection and all groups
        #
        self.update_all_groups()
        self.selected_groups=self.calcs[new_calc]['selected_groups']
        self.update_selected_groups()
        #
        # Set pdbfilename
        #
        self.pdbfilename=new_calc
        #
        # Read the PDB file
        #
        import Protool
        self.Protool_instance=Protool.structureIO()
        self.Protool_instance.parsepdb(lines=self.calcs[new_calc]['pdblines'])
        return

    #
    # -----
    #

    def update_display_calc(self,event=None):
        """Simply update the drawing of the titration curves"""
        self.update_selected_groups()
        return

    #
    # --------
    #

    def select_all(self,event=None):
        """Select all groups in the window"""
        count=0
        for group in self.groups:
            self.group_disp.selection_set(count)
            count=count+1
        return
    #
    # --------
    #

    def select_group(self):
        """Transfer the selected group to the 'selected' window"""
        for sel in self.group_disp.curselection():
            group_num=int(sel)
            sel_group=self.groups[group_num]
            if not sel_group in self.selected_groups:
                self.selected_groups.append(sel_group)
        self.update_selected_groups()
        return

    #
    # -----
    #
    
    def select_irregular_titcurvs(self):
        """Select irregular titration curves based on a fit to the Henderson-Hasselbalch equation"""
        self.selected_groups=[]
        X=pKanalyse()
        for group in self.groups:
            #
            # Include only groups with a pKa value
            #
            if self.pkavals[group]['pKa']:
                solution,sq=X.fit_to_henderson(self.titration_curves[group])
                if abs(solution[0])<self.irr_cutoff.get():
                    if not group in self.selected_groups:
                        self.selected_groups.append(group)
        self.update_selected_groups()
        return

    #
    # --------
    #

    def display_details_sel(self,event=None):
        """Display the details of the groups clicked in the selected window"""
        self.details.delete(1.0, END)
        #
        # Loop over all groups
        #
        for sel in self.selected.curselection():
            group_num=int(sel)
            group=self.selected_groups[group_num]
            self.add_group_to_details(group)
        return

    #
    # --------
    #

    def display_details(self,event=None):
        """Display the details of the groups clicked in the list of all groups"""
        #
        # Do the details
        #
        self.details.delete(1.0, END)
        #
        # Loop over all groups
        #
        for sel in self.group_disp.curselection():
            group_num=int(sel)
            group=self.groups[group_num]
            self.add_group_to_details(group)
        return

    #
    # -------
    #


    def add_group_to_details(self,group):
        """Add the details of a single group to the details window"""
        X=pKarun.pKa_general.pKanalyse()
        order=['modelpK','desolv','backgr','intpKa','delec','dpKa','pKa']
        unit_dict={'modelpK':'','desolv':'pKa units','backgr':'pKa units','delec':'pKa units','dpKa':'pKa units','pKa':'','intpKa':''}
        self.details.insert(END,'%-10s: %s\n' %('Name',group))
        for key in order:
            unit=unit_dict[key]
            if key=='intpKa':
                try:
                    self.pkavals[group]['intpKa']=self.pkavals[group]['modelpK']+self.pkavals[group]['desolv']+self.pkavals[group]['backgr']
                except:
                    self.pkavals[group]['intpKa']=None
            elif key=='dpKa':
                if self.pkavals[group]['pKa']:
                    self.pkavals[group]['dpKa']=self.pkavals[group]['pKa']-self.pkavals[group]['modelpK']
                else:
                    self.pkavals[group]['dpKa']=None
            #
            # Insert the text
            #
            if self.pkavals[group][key]:
                text='%-10s: %5.1f %s\n' %(key,self.pkavals[group][key],unit)
            else:
                text='%-10s: %5s %s\n' %(key,'Unknown',unit)
            self.details.insert(END,text)
        #
        # Insert the fit to the Henderson-Hasselbalch equation
        #
        solution,sq=X.fit_to_henderson(self.titration_curves[group])
        self.details.insert(END,'Fit to Henderson-Hasselbalch\na %.2f (1 for perfect HH)\npKa value from HH fit: %5.2f\n' %(abs(solution[0]),abs(solution[1])))
        self.details.insert(END,'---------------------\n\n')
        return

    #
    # -----
    #

    def show_details_table(self,event=None):
        """Open a new window with all details for the groups selected"""
        details_win=Toplevel()
        X=pKanalyse()
        row=0
        order=['Group','modelpK','desolv','backgr','intpKa','delec','dpKa','pKa','HH fit','HH fit pKa']
        #
        # Insert table headings
        #
        column=0
        for key in order:
            tl=Label(details_win,text=key,anchor='e',borderwidth=2,width=10,relief='ridge')
            tl.grid(row=row,column=column)
            column=column+1
        row=row+1
        #
        # Insert the data
        #
        for group in self.selected_groups:
            unit_dict={'modelpK':'','desolv':'pKa units','backgr':'pKa units','delec':'pKa units','dpKa':'pKa units','pKa':'','intpKa':''}
            tl=Label(details_win,text=group,anchor='e',borderwidth=2,relief='ridge',width=10)
            tl.grid(row=row,column=0)
            column=1
            for key in order:
                #
                # We've already inserted the name, and we will insert the HH parms later
                #
                if key in ['Group','HH fit','HH fit pKa']:
                    continue
                #
                # Inser the rest
                unit=unit_dict[key]
                if key=='intpKa':
                    try:
                        self.pkavals[group]['intpKa']=self.pkavals[group]['modelpK']+self.pkavals[group]['desolv']+self.pkavals[group]['backgr']
                    except:
                        self.pkavals[group]['intpKa']=None
                elif key=='dpKa':
                    if self.pkavals[group]['pKa']:
                        self.pkavals[group]['dpKa']=self.pkavals[group]['pKa']-self.pkavals[group]['modelpK']
                    else:
                        self.pkavals[group]['dpKa']=None
                #
                # Insert the text
                #
                if self.pkavals[group][key]:
                    text='%5.1f' %(self.pkavals[group][key])
                else:
                    text='%5s' %('Unknown')
                #
                # Fill the cell
                #
                tl=Label(details_win,text=text,anchor='e',borderwidth=2,width=10,relief='ridge')
                tl.grid(row=row,column=column)
                column=column+1
            #
            # Insert the fit to the Henderson-Hasselbalch equation
            #
            solution,sq=X.fit_to_henderson(self.titration_curves[group])
            tl=Label(details_win,text='%5.2f' %(abs(solution[0])),anchor='e',borderwidth=2,width=10,relief='ridge')
            tl.grid(row=row,column=column)
            column=column+1

            tl=Label(details_win,text='%5.2f' %(abs(solution[1])),anchor='e',borderwidth=2,width=10,relief='ridge')
            tl.grid(row=row,column=column)
            #
            row=row+1
        
        return

    #
    # ----
    #
    
    def find_groups(self):
        #
        # Get perturbed pKas
        #
        self.selected_groups=[]
        type=self.type.get()
        qualifier=self.qualifier.get()
        amount=self.amount.get()
        #
        # Find the groups
        #
        import pKarun
        for group in self.groups:
            pka=self.pkavals[group]['pKa']
            model=self.pkavals[group]['modelpK']
            this_type=isacid(group)
            resname=get_resname(group)
            #
            # Exclude types
            #
            if type=='Acids' and this_type!=1:
                continue
            elif type=='Bases' and this_type==1:
                continue
            elif type=='Asp + Glu' and resname!='GLU' and resname!='ASP':
                continue
            elif type=='Lys + Arg' and resname!='ARG' and resname!='LYS':
                continue
            elif type in ['HIS','ARG','CYS','ASP','GLU','TYR','LYS']:
                if resname!=type:
                    continue
            #
            # Get the right qualifier
            #
            if pka and model:
                if self.criteria_fulfilled(qualifier,amount,pka,model):
                    if not group in self.selected_groups:
                        self.selected_groups.append(group)
                else:
                    # Not added
                    continue
            else:
                if not group in self.selected_groups:
                    self.selected_groups.append(group)
        self.update_selected_groups()
        return

    #
    # ----
    #

    def criteria_fulfilled(self,qualifier,value,test_value,test_value_reference=None):
        #
        # Tests if a certain criteria is fulfilled
        #
        if qualifier=='equal to':
            if value==test_value:
                return 1
            return None
        elif qualifier=='raised by more than':
            if test_value-test_value_reference>value:
                return 1
            return None
        elif qualifier=='lowered by more than':
            if test_value-test_value_reference<value:
                return 1
            return None
        elif qualifier=='greater than':
            if test_value>value:
                return 1
            return None
        elif qualifier=='smaller than':
            if test_value<value:
                return 1
            return None
        else:
            raise 'Unknown qualifier',qualifier

    #
    # -------
    #

    def update_selected_groups(self):
        """Update the display in the main window - both the list of groups, and the drawing of titration curves"""
        #
        # Sort
        #
        self.selected_groups.sort()
        self.selected.delete(0, END)
        for group in self.selected_groups:
            self.selected.insert(END,group)
        #
        # Copy it to the master array
        #
        self.calcs[self.current_selected_calc]['selected_groups']=self.selected_groups[:]
        #
        # Redraw the titration curves
        #
        curves={}
        for calc in self.calcs.keys():
            if self.calc_disp[calc].get()==1:
                tit_curvs=self.calcs[calc]['titcurv']
                for group in self.calcs[calc]['selected_groups']:
                    label=self.calcs[calc]['label']
                    try:
                        curves[label+':'+group]=tit_curvs[group].copy()
                    except:
                        print 'Could not find group'
                        print self.titration_curves.keys()
                        raise Exception('Could not find group!: %s' %group)
        #
        # Update the drawing of curves
        #
        self.update_curves(curves,labels=1)
        #
        # If we are showing the stability curve, then update that one too
        #
        stab_status=self.stability_var.get()
        if stab_status=='on':
            self.do_stab_curve()
        return
    #
    # ----
    #

    def remove_group(self):
        #
        # Remove a group from the selected list
        #
        remove=[]
        for sel in self.selected.curselection():
            group_num=int(sel)
            remove.append(group_num)
        remove.sort()
        remove.reverse()
        for group_num in remove:
            self.selected_groups=self.selected_groups[:group_num]+self.selected_groups[group_num+1:]
        self.update_selected_groups()
        return   
    #
    # ---------
    #

    def show_system(self,other_groups_model=None,pH=7.0,pKasys_mode='Boltzmann'):
        """Show the subsystem in pKa_system"""
        #
        # If no groups selected then do nothing
        #
        if len(self.selected_groups)==0:
            return
        #
        # If we have too many groups then complain
        #
        if len(self.selected_groups)>10 and pKasys_mode!='Boltzmann (C++)':
            import tkMessageBox
            tkMessageBox.showwarning(
                "Too many groups",
                "You have selected too many groups for interactive analysis (%d)" %(len(self.selected_groups)))
            return
        #
        # Open dialog to let the user choose how to model the rest of the protein
        #
        self.pH_entry=StringVar()
        self.pH_entry.set(str(pH))
        self.choice=IntVar()
        if len(self.selected_groups)<len(self.pkavals.keys()):
            if other_groups_model is None:
                choice_win=Toplevel()
                self.set_geometry(self.master,choice_win)
                choice_win.title('Choose model for rest of protein')
                Label(choice_win,
                      text='How do you want to model the effect of all other groups in pKa_system?').grid(row=0,column=0,columnspan=2)
                self.choice.set(0)
                Radiobutton(choice_win,
                            text='Ignore all other groups',
                            value=0,
                            variable=self.choice).grid(row=1,column=0,sticky='news',columnspan=2)
                Radiobutton(choice_win,
                            text='Model as non-titrating point charges using charge at pH',
                            value=1,variable=self.choice).grid(row=2,column=0,sticky='news',columnspan=1)
                
                Entry(choice_win,textvariable=self.pH_entry,width=5,bg='white').grid(row=2,column=1,sticky='news')
                Radiobutton(choice_win,
                            text='Model as fixed titrating groups',
                            value=2,
                            variable=self.choice).grid(row=3,column=0,sticky='news',columnspan=2)
                Button(choice_win,text='Continue',command=choice_win.destroy,fg='green').grid(row=4,column=0,sticky='news')
                cancel=None
                def cancel_command():
                    #
                    # Function to be called if teh user cancels
                    #
                    cancel=1
                    choice_win.destroy()
                    return
                Button(choice_win,text='Cancel',command=cancel_command,fg='red').grid(row=4,column=1,sticky='news')
                self.master.wait_window(choice_win)
                if cancel:
                    return
            else:
                #
                # We have input from the function call
                #
                self.choice.set(other_groups_model)
            non_system_model=self.model_rest_of_protein(self.choice.get(),self.pH_entry.get())
        else:
            non_system_model=None
        #
        # Open pKa_system with these groups
        #
        self.matrix=self.X.read_matrix()
        S=pKa_system.pKa_system(numgroups=len(self.selected_groups),parent_application=self.master,update=False)
        S.non_system_groups=non_system_model
        self.pKasys=S
        self.pKasys.pkamethod_sel.set(pKasys_mode)
        self.pKasys.update_pkasystem_curves()
        #
        # Set all the interactions
        #
        id=0
        #
        # Tell pKa system to not update the view while we set parameters
        #
        S.stab_test_on=1
        #
        for group in self.selected_groups:
            #
            # Set the acid/base value
            #
            import pKarun
            if isacid(group):
                S.groups[id].acid_base.set(1)
            else:
                S.groups[id].acid_base.set(0)
            #
            # Set the name
            #
            S.groups[id].name.set(group)
            #
            # Set the model pKa value
            #
            this_grp=self.pkavals[group]
            S.groups[id].modelpK=this_grp['modelpK']
            #
            # Get the intrinsic pKa for this one
            #
            intpka=this_grp['modelpK']+this_grp['backgr']+this_grp['desolv']
            S.groups[id].intpka.set(intpka)
            S.groups[id].intpka_scale.set(intpka)
            #
            # Now for the interaction energies
            #
            for group2_id in S.groups[id].intenes.keys():
                group1_id=id
                group1_name=self.selected_groups[group1_id]
                group2_name=self.selected_groups[group2_id]
                #
                # Get the interaction energy
                #
                record=self.matrix[group1_name][group2_name]
                intene=abs(record[0]-record[1]-record[2]+record[3])
                #
                # Set it - we always use absolute values in pKa_system
                #
                S.groups[id].intenes_scales[group2_id].set(intene)         
            #
            # Update the id
            #
            id=id+1
        #
        # Return pKa_system to normal mode
        #
        self.master.update()
        S.stab_test_on=None
        S.update_pkasystem_curves()
        return
    
    #
    # ----
    #

    def model_rest_of_protein(self,model,set_pH):
        """Compile the information for modelling the rest of teh protein
        model 0 = ignore
        model 1 = non-titrating point charges (select charge at pH 7.0)
        model 2 = fixed titrating charges
        """
        wait_win=Toplevel()
        self.set_geometry(self.master,wait_win)
        wait_win.title('Please wait')
        Label(wait_win,text='Calculating cumulated effect of all other groups').grid(row=0,column=0)
        self.master.update_idletasks()
        wait_win.update_idletasks()
        if model==0:
            wait_win.destroy()
            return None
        elif model==2 or model==1:
            #
            # Initialise the array that holds information on all other groups
            #
            non_system_groups={}
            #
            # Get the interactions with all other groups at each pH value
            #
            group_id=0
            for group_name in self.selected_groups:
                #
                #
                non_system_groups[group_id]={}
                for pH in range(0,201):
                    real_pH=float(pH)/10.0
                    sum_interaction=0.0
                    for group2_name in self.matrix[group_name].keys():
                        if group2_name in self.selected_groups:
                            #
                            # Don't include the selected groups
                            #
                            pass
                        else:
                            #
                            # For all other groups calculate the interaction at each pH
                            #
                            record=self.matrix[group_name][group2_name]
                            intene=record[0]-record[1]-record[2]+record[3]
                            #
                            # Get the charge for this group at this pH
                            #
                            if model==2:
                                #
                                # Get charge at the right pH
                                #
                                charge=self.get_charge(group2_name,real_pH)
                            elif model==1:
                                #
                                # Just grab charge at the pH specified
                                #
                                charge=self.get_charge(group2_name,float(set_pH))
                            else:
                                print 'Model undefined'
                                raise Exception
                            if charge is None:
                                import tkMessageBox
                                tkMessageBox.showwarning("Could not find charge",
                                                         "Could not find charge for %s at pH %5.2f\nDid you calculate charges for a wide enough pH range?" %(group2_name,real_pH),
                                                         parent=self.master)
                            sum_interaction=sum_interaction+intene*abs(charge) # The intene has the right sign already
                    #
                    # Store the info
                    #
                    non_system_groups[group_id][real_pH]=sum_interaction
                #
                # Update id
                #
                group_id=group_id+1
        wait_win.destroy()
        return non_system_groups

    #
    # ---------
    #

    def goto_pKa_system(self):
        """Open dialog asking for how many groups"""
        import tkSimpleDialog
        num_groups=tkSimpleDialog.askinteger("Number of groups","Please enter number of titratable groups",
                                             parent=self.master)
        S=pKa_system.pKa_system(num_groups,self.master)
        return

    #
    # --------
    #

    def select_strongly_interacting_groups(self):
        """Select groups that interact strongly with the selected groups"""
        self.matrix=self.X.read_matrix()
        append_groups={}
        #
        # Loop over all selected groups and check for interaction energies
        #
        interaction_energy_cutoff=self.int_cutoff.get()
        for group_sel in self.selected_groups:
            #self.selected_groups.append(group)
            #
            # Loop over all groups
            #
            for group in self.groups:
                record=self.matrix[group_sel][group]
                intene=record[0]-record[1]-record[2]+record[3]
                if abs(intene)>=interaction_energy_cutoff:
                    append_groups[group]=1
        #
        # Found all, now add them
        #
        for group in append_groups.keys():
            if not group in self.selected_groups:
                self.selected_groups.append(group)
        self.update_selected_groups() 
        return

    #
    # ------------
    #


    def pka_calc(self):
        self.master.destroy()
        X=Main_pKa_Calc_Window().mainloop()
        return

    #
    # -----
    #
    
    def pKa_system(self):
        self.master.destroy()
        X=pKa_system.pKa_system().mainloop()
        return


    #
    # -------
    #

    def load_titration(self,pdbfilename=None):
        """Load the results of a pKa calculation"""
        import os
        if not pdbfilename:
            pdbfilename=tkFileDialog.askopenfilename(defaultextension='.pdb',
                                                     initialdir=os.getcwd(),
                                                     filetypes=[("WHAT IF pKa dir","*.pdb"),
                                                                ("pKaTool pKa calc file","*.ppc"),
                                                                ("MEAD pKa file","*.pka"),
                                                                ("All files","*.*")],
                                                     parent=self.master)
        #
        # Did we get a name?
        #
        if not pdbfilename:
            return
        pdbfilename=os.path.join(os.getcwd(),pdbfilename)
        if not os.path.isfile(pdbfilename):
            tkMessageBox.showwarning("File not found",
                                    'File not found: %s' %pdbfilename
                                    ,parent=self.master)
            return
        #
        # Is it a WHAT IF type file?
        #
        if pdbfilename[-4:]=='.pdb':
            #
            # Open the file and assess status
            #
            import os
            pdbfilename_short=os.path.split(pdbfilename)[1]
            label=self.labels[self.label_count]
            self.calcs[pdbfilename_short]={'instance':pKaIO.pKaIO(pdbfilename),'selected_groups':[],'label':label}
            self.X2=self.calcs[pdbfilename_short]['instance']
            if not self.X2.assess_status():
                tkMessageBox.showwarning("pKa calculations not complete",
                                         'Please perform a complete pKa calculation for this PDB file'
                                         ,parent=self.master)
                del self.calcs[pdbfilename_short]
                return
            #
            # Read the titration curve
            #
            self.calcs[pdbfilename_short]['titcurv']=self.X2.readtitcurv()
            #
            # read the pKa values
            #
            self.calcs[pdbfilename_short]['pkavals']=self.X2.readpka()
            #
            # Store the PDB file
            #
            fd=open(pdbfilename)
            self.calcs[pdbfilename_short]['pdblines']=fd.readlines()
            fd.close()
            #
            # Add to the display pulldown menu
            #
            self.label_count=self.label_count+1
            self.add_to_pulldown(pdbfilename_short,label)
            self.calc_select.set(pdbfilename_short)
            self.calc_disp[pdbfilename_short].set(1)
            self.update_selected_calc()
        elif pdbfilename[-4:]=='.ppc':
            #
            # pKaTool file
            #
            prog=Toplevel()
            prog.geometry('250x100')
            prog.title('Load progress')
            self.set_position(prog)
            self.load_file=StringVar()
            self.load_file.set('Opening file...')
            Label(prog,textvariable=self.load_file).grid(row=0,column=0)
            prog.update_idletasks()


            import pickle, copy, os
            fd=open(pdbfilename)
            calcs=pickle.load(fd)
            fd.close()

            import os
            count=0
            for calc in calcs.keys():
                count=count+1
                self.load_file.set('Parsing %s (%d of %d)' %(os.path.split(calc)[1],count,len(calcs.keys())))
                prog.update_idletasks()
                calc_name=os.path.split(calc)[1]
                self.calcs[calc_name]=copy.deepcopy(calcs[calc])
                self.label_count=self.label_count+1
                self.add_to_pulldown(calc_name,'')
                self.calc_select.set(calc_name)
                self.calc_disp[calc_name].set(1)
                self.update_selected_calc()
            prog.destroy()
        elif pdbfilename[-4:]=='.pka':
            #
            # MEAD file
            #
            import os
            pdbfilename_short=os.path.split(pdbfilename)[1]
            label=self.labels[self.label_count]
            self.calcs[pdbfilename_short]={'instance':pKaIO.pKaIO(pdbfilename),'selected_groups':[],'label':label}
            self.X2=self.calcs[pdbfilename_short]['instance']
            #
            # Read the titration curve
            #
            self.calcs[pdbfilename_short]['titcurv']=self.X2.readtitcurv(pdbfilename)
            #
            # read the pKa values
            #
            self.calcs[pdbfilename_short]['pkavals']=self.X2.readpka(pdbfilename)
            #
            # Put the pKa value in the titration curves
            #
            for group in self.calcs[pdbfilename_short]['titcurv']:
                self.calcs[pdbfilename_short]['titcurv'][group]['pKa']=self.calcs[pdbfilename_short]['pkavals'][group]['pKa']
            #
            # Store the PDB file
            #
            fd=open(pdbfilename)
            lines=fd.readlines()
            fd.close()
            pdblines=self.X2.read_MEADfile(lines,mode='pdb')
            self.calcs[pdbfilename_short]['pdblines']=pdblines
            #
            # Add to the display pulldown menu
            #
            self.label_count=self.label_count+1
            self.add_to_pulldown(pdbfilename_short,label)
            self.calc_select.set(pdbfilename_short)
            self.calc_disp[pdbfilename_short].set(1)
            self.update_selected_calc()
        return

    #
    # -----
    #

    def save_titration(self):
        """Save all titrations in a pickled file"""
        import tkFileDialog, os, pickle
        savedir=os.getcwd()
        filename=tkFileDialog.asksaveasfilename(defaultextension=".ppc",
                                                initialfile='calculation1.ppc',
                                                initialdir=savedir,
                                                filetypes=[("pKaTool pKa calc file","*.ppc"),
                                                           ("WHAT IF pKa dir format","*.pdb"),
                                                           ("All files","*.*")],
                                                parent=self.master)
        if not filename:
            return
        if filename[-4:]=='.ppc':
            fd=open(filename,'w')
            pickle.dump(self.calcs,fd)
            fd.close()
        elif filename[-4:]=='.pdb':
            print 'Not implemented yet'
            return
        else:
            #
            # Save in ppc format anyway
            #
            fd=open(filename,'w')
            pickle.dump(self.calcs,fd)
            fd.close()
        return

    #
    # -----
    #

    def subtract_titration(self):
        """Subtract the second pKa calculation from the first to find differences in titration curves"""
        #
        # Choose the calc to subtract from
        #
        self.choose_win=Toplevel()
        lbl=Label(self.choose_win,text='Which calculation should we subtract from %s' %self.calc_select.get())
        lbl.grid(row=0,column=0,columnspan=4)
        #
        self.sub_calc=StringVar()
        self.sub_calc.set('Choose a dataset')
        s_button=Menubutton(self.choose_win,textvariable=self.sub_calc,relief=RAISED)
        # the menu
        s_menu=Menu(s_button,tearoff=0)
        s_button['menu']=s_menu
        for name in self.calcs.keys():
            s_menu.add_radiobutton(label=name,variable=self.sub_calc,value=name,indicatoron=1)
        s_button.grid(row=1,column=0,sticky='news')
        #
        Button(self.choose_win,text='Subtract',command=self.destroy_choose).grid(row=2,column=0)
        Button(self.choose_win,text='Cancel',command=self.cancel_subtract).grid(row=2,column=2)
        self.master.wait_window(self.choose_win)
        #
        # Create new set
        #
        first=self.current_selected_calc
        second=self.sub_calc.get()
        name=self.calcs[first]['label']+'-'+self.calcs[second]['label']
        #
        # Subtract the pKa values
        #
        intpkas=[]
        sub_pka={}
        for group in self.calcs[first]['pkavals'].keys():
            if self.calcs[second]['pkavals'].has_key(group):
                sub_pka[group]={}
                for property in self.calcs[first]['pkavals'][group].keys():
                    if self.calcs[second]['pkavals'][group].has_key(property):
                        diff=self.calcs[first]['pkavals'][group][property]-self.calcs[second]['pkavals'][group][property]
                        sub_pka[group][property]=diff
        #
        # Subtract the titration curves
        #
        sub_tits={}
        for group in self.calcs[first]['titcurv'].keys():
            if self.calcs[second]['titcurv'].has_key(group):
                sub_tits[group]={}
                for pH in self.calcs[first]['titcurv'][group].keys():
                    if self.calcs[second]['titcurv'][group].has_key(pH):
                        diff=self.calcs[first]['titcurv'][group][pH]-self.calcs[second]['titcurv'][group][pH]
                        sub_tits[group][pH]=diff
        #
        # Done subtracting. Create new set
        #
        label=self.labels[self.label_count]
        self.label_count=self.label_count+1
        self.calcs[name]={'titcurv':sub_tits.copy(),'label':label,'instance':None,'selected_groups':[],'pkavals':sub_pka.copy()}
        #
        # Add to the display pulldown menu
        #
        self.add_to_pulldown(name,label)
        self.calc_select.set(name)
        self.calc_disp[name].set(1)
        self.update_selected_calc()
        return

    #
    # -----
    #

    def destroy_choose(self,event=None):
        """Destroy the choose_window"""
        self.choose_win.destroy()
        return

    #
    # -----
    #

    def cancel_subtract(self,event=None):
        """Cancel the subtraction action"""
        self.sub_calc.set('')
        self.destroy_choose()
        return


    #
    # ----
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
    # ----
    #
        
    def set_geometry(self,pwidget,widget):
        """Set the position of widget in the middle of pwidget"""
        w,h,x,y=self.get_geometry(pwidget)
        pwidget.update()
        widget.update()
        sw,sh,dummy,dummy2=self.get_geometry(widget)
        xoffset=int((w-sw)/2)
        yoffset=int((h-sh)/2)
        widget.geometry('+%d+%d' %(x+xoffset,y+yoffset))
        return

#
# -----
#



def main():
    """pKaTool main function"""
    import sys
    X=Main()
    #
    # If a filename was given then load that file
    #
    if len(sys.argv)==2:
        filename=sys.argv[1]
        import os
        if os.path.isfile(filename):
            X.load_titration(filename)
    #
    # Enter Tk mainloop
    #
    X.mainloop()
    return

if __name__=="__main__":
    main()

    
