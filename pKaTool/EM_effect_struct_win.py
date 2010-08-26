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

from Tkinter import *
import pKa_base
import EM_effect

class EM_effect_struct(Frame,pKa_base.pKa_base,EM_effect.EM_effect):

    def __init__(self,parent=None):
        """Init the class"""
        self.parent_application=None
        if parent:
            self.parent_application=parent
        #
        # Open the window
        #
        if not self.parent_application:
            Frame.__init__(self)
            self.EM_window=self.master
        else:
            self.EM_window=Toplevel()
        self.EM_window.geometry('+200+600')
        self.EM_window.title('Structure/sequence display')
        #
        # Draw the pH vs. chem shift graph
        #
        self.linewidth=2
        self.titwidth=1200
        self.titheight=450
        self.tc=Canvas(self.EM_window,bd=5,bg='white',width=self.titwidth,
                       height=self.titheight,
                       scrollregion=(0,0,self.titwidth,self.titheight))
        self.tc.xview("moveto", 0)
        self.tc.yview("moveto", 0)
        self.tc.grid(row=0,column=0)
        #
        # Axes
        #
        self.set_graph_size(x_max=129.0,x_min=0.0,y_max=12.0,y_min=0.0,
                            x_tick_level=20,y_tick_level=2.0,width=self.titwidth,height=self.titheight)
        self.draw_ordinates(self.tc,y_label='pH',x_label='Residue')
        #
        # Pulldown menu
        #
        self.menu=Menu(self.EM_window)
        #
        # File menu
        #
        self.file_menu=Menu(self.menu,tearoff=0)
        self.file_menu.add_command(label='Load PDB file',command=self.load_pdb)
        self.file_menu.add_command(label='Load Ekin titration data',command=self.load_full_titdata)
        self.file_menu.add_command(label='Exit',command=self.quit)
        self.menu.add_cascade(label='File',menu=self.file_menu)
        #
        # Configure the menu
        #
        self.EM_window.config(menu=self.menu)
        return
    #
    # ----
    #

    def draw_titrations(self,titdata):
        """Draw all the titrations"""
        residues=titdata.keys()
        residues.sort()
        for residue in residues:
            resnumber=int(residue[1:])
            if len(titdata[residue])==0:
                continue
            last_pH,last_chemshift=titdata[residue][0]
            for pH,chemshift in titdata[residue][1:]:
                x,y=self.get_xy(resnumber,last_pH)
                x1,y1=self.get_xy(resnumber+1,pH)
                #
                # Get the colour - depends on the chemical shift change
                #
                colour=self.get_colour_from_chemshift((chemshift-last_chemshift)/(pH-last_pH))
                self.tc.create_rectangle(x,y,x1,y1,fill=colour,outline=colour)
                #
                # Update last pH and chemshift
                #
                last_chemshift=chemshift
                last_pH=pH
        return

    #
    # -----
    #

    def get_colour_from_chemshift(self,dchemshift):
        """Get the colour corresponding to a certain chemshift change"""
        import numpy, string
        col_min=numpy.array([100,100,100])
        if dchemshift>0.0:
            col_max=numpy.array([0,0,255])
        elif dchemshift<=0.0:
            col_max=numpy.array([255,0,0])
        #
        # Calculate the colour
        #
        dchemshift=abs(dchemshift)
        max_change=0.5
        if dchemshift>max_change:
            dchemshift=max_change
        colour=(col_max-col_min)/max_change*dchemshift+col_min
        text='#'
        for col in colour:
            text=text+string.zfill(hex(int(col))[2:],2)
        return text
        
        
if __name__=='__main__':
    X=EM_effect_struct()
    #X.load_pdb('2lzt.pka.pdb')
    X.mainloop()
