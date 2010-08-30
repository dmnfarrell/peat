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

from Tkinter import *
import tkFileDialog
import Pmw

class OptparseDialog(object):
    """Displays an automatically generated dialog using
       an optparse object"""

    def __init(self):
        return

    def createWidgets(self, parent, parser):
        """Create widgets for parser options"""
        #CheckBox for boolean options
        #ComboBox for "choice" options
        self.frame=Frame(parent)
        bal=Pmw.Balloon(self.frame)
        self.optvars={}
        optvars = self.optvars
        r=0
        for option in parser.option_list:
            if option.dest is None:
                continue
            if option.help is None:
                option.help = u''
            #print option.action, option.type,  option.default,option.metavar
            l=Label(self.frame, text=option.dest)
            l.grid(row=r,column=0, sticky='e')
            bal.bind(l, option.help)
            if 'store' == option.action:
                if option.type == 'string':
                    optvars[option.dest] = StringVar()
                    Entry(self.frame, textvariable=optvars[option.dest]).grid(row=r,column=1)

                if option.metavar == 'FILE':
                    def func(var):
                        def getfile():
                            fn=tkFileDialog.askopenfilename(parent=self.frame)
                            var.set(fn)
                        return getfile
                    Button(self.frame,text='..',
                            command=func(optvars[option.dest])).grid(row=r,column=2)

            elif option.action in ('store_true', 'store_false'):
                optvars[option.dest]=BooleanVar()
                optvars[option.dest].set(option.default)
                Checkbutton(self.frame, variable=optvars[option.dest]).grid(row=r,column=1)

            r+=1
        return self.frame, optvars
