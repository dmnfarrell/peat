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

class group_control:
    """Defines a single group and the interaction energies
    with all other groups"""
    
    def __init__(self,parent,row_num,id_num,total,int_ene=1,colour_order=['black'],window=None):
        """Initialise all parameters and arrays"""
        #
        # Windows
        #
        self.parent=parent
        self.win=self.parent.win
        if window:
            self.win=window
        #
        # pKa related arrays
        #
        self.modelpK=None
        self.intenes_required=int_ene
        #
        # Keep track of colours
        #
        self.colour_order=colour_order[:]
        self.colour=colour_order[id_num%len(colour_order)]
        self.fg_col='black'
        #
        # Should we do the labels?
        #
        if id_num==0:
            row_num=self.labels(total,row_num)
        else:
            row_num=row_num+1
        #
        # Draw a single set of controls for one group
        #
        self.column=0
        #
        # Define all the variables
        #
        self.active=IntVar()
        self.acid_base=IntVar()
        self.acid_base_text=StringVar()
        self.intpka=DoubleVar()
        self.name=StringVar()
        self.style=IntVar()
        self.style_text=StringVar()
        #
        # Active?
        #
        self.active_btn=Checkbutton(self.win,
                                    variable=self.active,  
                                    bg=self.colour,
                                    fg=self.fg_col,
                                    activebackground=self.colour,
                                    selectcolor='white',
                                    disabledforeground='white',
                                    indicatoron=1)
        self.active_btn.grid(row=row_num,
                             column=self.column,
                             sticky='news')
        self.column=self.column+1
        #
        # Name
        #
        name=Label(self.win,textvariable=self.name,
                   background=self.colour,
                   foreground=self.fg_col,
                   activebackground=self.colour,
                   activeforeground=self.fg_col)
        name.grid(row=row_num,column=self.column,sticky='news')
        name.bind('<ButtonRelease-1>',self.change_name)
        self.column=self.column+1
        #
        # Label for group number 
        #
        self.no_label=Label(self.win,
                            text='%d' %id_num,
                            bg=self.colour,
                            fg=self.fg_col)
        self.no_label.grid(row=row_num,
                           column=self.column,
                           sticky='news')
        #
        # Acid/base button
        #
        self.column=self.column+1
        self.ab_text_column=self.column
        self.ab_button=Checkbutton(self.win,
                                   textvariable=self.acid_base_text,
                                   variable=self.acid_base,
                                   bg=self.colour,
                                   fg=self.fg_col,
                                   activebackground=self.colour,
                                   activeforeground=self.fg_col)
        self.ab_button.grid(row=row_num,column=self.column,sticky='news')
        #
        # Intrinsic pKa
        #
        self.column=self.column+1
        self.intpka_scale=Scale(self.win,from_=0.0,to=14.0,resolution=0.1,
                                orient='horizontal',
                                relief='ridge',
                                bg=self.colour,
                                activebackground=self.colour,
                                fg=self.fg_col)
                  
        self.intpka_scale.grid(row=row_num,column=self.column)
        #
        # Binding for fixing this energy during minimisation
        #
        self.intpka_scale.bind('<Button-3>',self.toggle_fix_intpka)
        #
        # Interaction energies
        #
        if self.intenes_required:
            self.intenes={}
            self.intenes_scales={}
            for num in range(id_num):
                self.column=self.column+1
                self.intenes_scales[num]=Scale(self.win,
                                               from_=0.0,to=20.0,
                                               resolution=0.1,
                                               orient='horizontal',
                                               relief='ridge')
                self.intenes_scales[num].grid(row=row_num,column=self.column)
                self.intenes[num]=DoubleVar()
                #
                # Binding for fixing energy
                #
                def toggle_fix_energy(event=None,scale=self.intenes_scales[num]):
                    """Toggle this energy fix"""
                    state=scale.cget('state')
                    if state=='disabled':
                        scale.configure(state='normal')
                        scale.configure(troughcolor='grey')
                    else:
                        scale.configure(state='disabled')
                        scale.configure(troughcolor='red')
                    return
                self.intenes_scales[num].bind('<Button-3>',toggle_fix_energy)
            #
            # Add '-'s for the interaction with self
            #
            self.column=self.column+1
            lbl=Label(self.win,text='-')
            lbl.grid(row=row_num,column=self.column)
            #
            # Empty cells for the rest
            #
            for num in range(total):
                if num<id_num+1:
                    continue
                #
                # Add an empty label 
                #
                self.column=self.column+1
                lbl=Label(self.win,text='')
                lbl.grid(row=row_num,column=self.column,sticky='news')
                self.win.columnconfigure(self.column,minsize=100)
            #
            # Add a field for the calculated pKa value
            #
            self.column=self.column+1
            self.pkavalue=StringVar()
            self.pka=Label(self.win,textvariable=self.pkavalue)
            self.pka.grid(row=row_num,column=self.column)
            self.win.columnconfigure(self.column,minsize=100)
            #
            # Label for quality of HH fit
            #
            self.column=self.column+1
            self.HHfit=StringVar()
            self.HH=Label(self.win,textvariable=self.HHfit)
            self.HH.grid(row=row_num,column=self.column)
            self.win.columnconfigure(self.column,minsize=50)
            #
            # Button for setting line style
            #
            self.column=self.column+1
            self.style_button=Checkbutton(self.win,
                                          textvariable=self.style_text,
                                          variable=self.style)
            self.style_button.grid(row=row_num,column=self.column)
        #
        # Set the values of all the variables
        #
        self.active_btn.select()
        self.acid_base.set(1)
        self.acid_base_text.set('Acid')
        self.intpka_scale.set(4.0)
        self.intpka.set(4.0)
        self.name.set('Group %d'%id_num)
        self.style_text.set('Line')
        self.style.set(1)
        #
        # done
        #
        return

    #
    # ----
    #

    def change_name(self,event=None):
        """Get new name"""
        import tkSimpleDialog
        name=tkSimpleDialog.askstring('Group name',prompt='Enter new name',
                                 initialvalue=self.name.get(),parent=self.win)
        if name:
            self.name.set(name)
        return

    #
    # -----
    #
    
    def toggle_fix_intpka(self,event=None):
        """Toggle the fix status of the intpka slider"""
        state=self.intpka_scale.cget('state')
        if state=='disabled':
            self.intpka_scale.configure(state='normal')
            self.intpka_scale.configure(troughcolor='grey')
        else:
            self.intpka_scale.configure(state='disabled')
            self.intpka_scale.configure(troughcolor='red')
        return

    #
    # ----------
    #

    def labels(self,total,row_num):
        """
        # Draw the labels
        """
        texts=['Active','Name','Group no.','Acid/Base','Int pKa']
        if self.intenes_required:
            for num in range(total):
                texts.append(str(num))
            texts.append('pKa 1/2')
            texts.append('HH fit pKa (slope/R^2)')
        texts.append('Style')
        #
        # Put the labels there
        #
        column_num=0
        for text_val in texts:
            if len(text_val)==1:
                colour=self.colour_order[int(text_val)%len(self.colour_order)]
            else:
                colour='#aaaaaa'
            fg_col='black'
            Label(self.win,text=text_val,bg=colour,fg=fg_col).grid(row=row_num,column=column_num,sticky='news')
            column_num=column_num+1
        return row_num+1

    #
    # ----------------
    #

    def update_scales(self):
        """Called from the parent class
        #
        # updates graph when scales are moved by user
        #
        """
        if not '%.1f'%self.intpka.get() == str(self.intpka_scale.get()):
            self.intpka.set(self.intpka_scale.get())
        #
        for num in self.intenes_scales:
            if not '%.1f'%self.intenes[num].get() == str(self.intenes_scales[num].get()):
                self.intenes[num].set(self.intenes_scales[num].get())
        import sys
        return

    #
    # ----
    #

    def update_scales_from_fit(self):
        """Called from the parent class
        #
        # updates scales when parameters are changed by fitter
        #
        """
        self.intpka_scale.set(self.intpka.get())
        for num in self.intenes_scales:
            self.intenes_scales[num].set(self.intenes[num].get())
        return

    #
    # -----
    #

    def update_group_control(self):
        """Update the acid/base label"""
        if self.acid_base.get()==1:
            self.acid_base_text.set('Acid')
        else:
            self.acid_base_text.set('Base')

        if self.style.get()==1:
            self.style_text.set('Line')
        else:
            self.style_text.set('Dashed')
        return
                

