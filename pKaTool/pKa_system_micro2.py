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

class Micro_states_helper:
    
    def toggle_kcat_window(self):
        """Open or close the kcat window"""
        if self.kcat_window_visible.get()==0:
            self.kcat_window.destroy()
        else:
            self.init_kcat_array()
            self.kcat_window=Toplevel()
            self.kcat_window.title('kcats for microstates')
            column=0
            for label in ['State','kcat for state']:
                Label(self.kcat_window,text=label).grid(row=0,column=column)
                column=column+1
            row=1
            for state in self.states:
                Label(self.kcat_window,text=state).grid(row=row,column=0)
                Entry(self.kcat_window,textvariable=self.kcat_microstates[state],width=5).grid(row=row,column=1)
                row=row+1
            #
            # Button for updating the plot
            #
            Button(self.kcat_window,command=self.update_pkasystem_curves,text='Update plot').grid(row=row+1,column=0)
        
        
        return
        
    #
    # ----
    #
    
    def init_kcat_array(self):
        """Initialise or complete the array holding kcat values for microstates"""
        self.init_states()
        if not hasattr(self,'kcat_microstates'):
            self.kcat_microstates={}
        for state in self.states:
            if not self.kcat_microstates.has_key(state):
                self.kcat_microstates[state]=DoubleVar()
                self.kcat_microstates[state].set(1.0)
        return
 
