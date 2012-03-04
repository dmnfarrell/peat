#!/usr/bin/env python
#
# DataPipeline - A data import and fitting tool
# Copyright (C) 2011 Damien Farrell
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
# Email: damien.farrell_at_ucd.ie
# Normal mail:
# Damien Farrell
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland
#

from Tkinter import *
import Images
import Pmw
from Base import Pipeline
from Editor import TextEditor

"""Implements a simple wizard in the application to allows users to
select and appropriate Importer and some other base settings"""

class HelperDialog(Frame):

    def __init__(self, parent=None):
        self.parent=parent
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
        else:
            self.main=Toplevel()
            self.master=self.main
        self.main.title('Configuration Helper')
        self.main.geometry('780x500+300+200')
        self.formatPanel()
        self.p = self.parent.p
        if self.parent:
            self.parent.wait_window(self.main)
        return

    def formatPanel(self):
        """Showing format options"""
        self.img1 = Images.formatsDiagram()
        label1 = Label(self.main,image=self.img1)
        label1.pack(fill=BOTH,padx=4,pady=4,ipadx=2,ipady=2)
        text = """The layout/format of your data is specified in the configuration
        file under the 'format' keyword. 8 typical formats
        are built into DataPipeline, covering most simple cases."""
        text=text.replace('\t','')
        Label(self.main,text=text,wraplength=600).pack(fill=Y,side=TOP,pady=4)
        panel2 = Frame(self.main)
        panel2.pack(fill=Y,side=BOTTOM)
        self.doButtons(panel2)
        return

    def doButtons(self, frame):
        """Buttons"""
        buttons=[['Show Config',self.showConfig],
                 ['Online Help',self.help],
                 ['Close',self.close]]

        for button,command in buttons:
            x=Button(frame,text=button,command=command)
            x.pack(side=LEFT,fill=BOTH,padx=2,pady=2)
        return

    def showConfig(self):
        """Show current config"""
        t = TextEditor(parent=self)
        print self.p.__dict__
        t.openFile(self.p.conffile)
        return

    def help(self):
        import webbrowser
        link='http://enzyme.ucd.ie/main/index.php/DataPipeline'
        webbrowser.open(link,autoraise=1)
        return

    def close(self,event=None):
        """Close the frame"""
        self.main.destroy()
        return
