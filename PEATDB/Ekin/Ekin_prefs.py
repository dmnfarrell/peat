#!/usr/bin/env python
#
# Protein Engineering Analysis Tool DataBase (PEATDB)
# Copyright (C) 2010 Damien Farrell & Jens Erik Nielsen
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

import sys,os
from Tkinter import *
import Pmw

class settings:
    """Class to provide dialog interface for ekin settings"""

    def __init__(self, parent=None):

        return

    #
    # Dialog to change ekin prefs
    #
    def show_settings(self):

        self.ekin_setupwin=Toplevel()
        self.ekin_setupwin.geometry('+300+450')
        self.ekin_setupwin.title('Ekin Settings')
        self.graphsizes = ['300x300', '500x400', '800x600', '1000x800']

        row=1
        lblfont=Label(self.ekin_setupwin,text='Font:')
        lblfont.grid(row=row,column=0,padx=3,pady=2)
        fontentry_button=Menubutton(self.ekin_setupwin,textvariable=self.font_input,
					relief=RAISED,width=16)
        fontentry_menu=Menu(fontentry_button,tearoff=0)
        fontentry_button['menu']=fontentry_menu
        #
        # Other fonts available
        #
        fts=['Arial','Courier','Verdana','Fixed','Times']
        for text in fts:
            #text='Font '+text
            fontentry_menu.add_radiobutton(label=text,
                                            variable=self.font_input,
                                            value=text,
                                            indicatoron=1)
        fontentry_button.grid(row=row,column=1, sticky='nes', padx=3,pady=2)

        row=3
        lblfontsize=Label(self.ekin_setupwin,text='Main Font Size:')
        lblfontsize.grid(row=row,column=0,padx=3,pady=2)
        fontsizeentry=Scale(self.ekin_setupwin,from_=8,to=20,resolution=1,orient='horizontal',
                            relief='ridge',variable=self.fontsize_input)

        fontsizeentry.grid(row=row,column=1, sticky='wens',padx=3,pady=2)
        row=4
        frame = Frame(self.ekin_setupwin)
        fontstyle_label = Label(frame, text='Font Style:')
        fontstyle_label.grid(row=0,column=0)
        fontstyle = Radiobutton(frame, text="plain", variable=self.fontstyle_input, value=0)
        fontstyle1 = Radiobutton(frame, text="bold", variable=self.fontstyle_input, value=1)
        fontstyle2 = Radiobutton(frame, text="italic", variable=self.fontstyle_input, value=2)
        fontstyle.grid(row=0,column=1)
        fontstyle1.grid(row=0,column=2)
        fontstyle2.grid(row=0,column=3)
        frame.grid(row=row,column=0,columnspan=2,sticky='news', padx=3,pady=2)
        row=5
        #listbox for canvas size
        self.graphsize_box = Pmw.ScrolledListBox(parent = self.ekin_setupwin,
                            labelpos='nw',
                            label_text='Graph Size:',
                            listbox_selectmode = 'single',
                            listbox_height = 4,
                            usehullsize = 1,
                            hull_width = 100,
                            hull_height = 100
        )
        self.graphsize_box.grid(row=row,column=0,sticky='nw')
        for item in self.graphsizes:
            self.graphsize_box.insert(END, item)
        self.graphsizes = map(int, self.graphsize_box.curselection())

        row=6
        c=Button(self.ekin_setupwin,text='Plotter settings', command=self.mainplotter.display_setup)
        c.grid(row=row,column=0,columnspan=2,sticky='news',padx=4,pady=4)

        row=7
        #
        # Apply Button
        #
        b = Button(self.ekin_setupwin, text="Apply Settings", command=self.apply_settings)
        b.grid(row=row,column=1,sticky='news',padx=4,pady=4)
        #
        # Close button
        #
        c=Button(self.ekin_setupwin,text='Close', command=self.close_ekin_setupwin)
        c.grid(row=row,column=0,sticky='news',padx=4,pady=4)
        #
        # Save Settings button
        #
        row=8
        c=Button(self.ekin_setupwin,text='Save as Default', command=self.save_preferences)
        c.grid(row=row,column=0,columnspan=2,sticky='news',padx=4,pady=4)
        return

    #
    # -----
    #
    def save_preferences(self):
        """Saves the sequence display settings to prefs file"""
        import preferences
        print 'Saving Ekin preferences'
        self.preferences.set('font',self.font_input.get())
        self.preferences.set('fontsize',self.fontsize_input.get())
        self.preferences.set('fontstyle',self.fontstyle_input.get())
        self.preferences.set('graphsize',self.graphsize)
        #self.preferences.set('showgraphlegend',self.showgraphlegend)
        #self.preferences.set('showgraphoptions',self.showgraphoptions)


        return
    #
    # -----
    #
    def load_preferences(self):
        """Loads the sequence display settings from prefs file, if present"""
        from PEATDB.Prefs import Preferences
        print 'Loading current Ekin preferences'
        self.preferences=Preferences('Ekin',{'canvas_height':600})

        try:
            f=self.preferences.get('font')
            self.font_input.set(f)
        except:
            self.preferences.set('font','fixed')
        try:
            f=self.preferences.get('fontsize')
            self.fontsize_input.set(f)
        except:
            self.preferences.set('fontsize',12)
        try:
            f=self.preferences.get('fontstyle')
            self.fontstyle_input.set(f)
        except:
            self.preferences.set('fontstyle',1)
        try:
            self.graphsize =self.preferences.get('graphsize')

        except:
            self.preferences.set('graphsize','500x400')
        return

    def close_ekin_setupwin(self):
        """Closes the display settings window"""
        self.ekin_setupwin.destroy()
        self.ekin_setupwin=None
        return

