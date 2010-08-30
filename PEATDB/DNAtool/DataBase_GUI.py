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

class Main(Frame):

    def __init__(self):
        #
        # Set up the main database window
        #
        # We need options for: Open database, new project
        #
        self.font='Times 12 bold'
        Frame.__init__(self)
        self.master.title("EATDB - Main Window")
        self.master.geometry("800x600")
        #
        # Welcome text
        #
        label1=Label(self.master, text="Welcome to EATDB",font=self.font)
        label1.grid(row=0,column=0, sticky=W)
        #
        # Add the pulldown menu
        #
        self.main_pulldown()
        
        self.mainloop()
        return


    def main_pulldown(self):
        #
        # Create the main pulldown menu
        #
        menu=Menu(self.master)
        #
        # File menu
        #
        filemenu=Menu(menu,tearoff=0)       
        filemenu.add_command(label='New Project',command=self.new_project)
        filemenu.add_command(label='Open',command=self.open_project)
        filemenu.add_separator()
        filemenu.add_command(label='Quit',command=self.quit)
        #
        # Add the file menu as a pulldown of the file menu
        #
        menu.add_cascade(label='File',menu=filemenu)
        self.master.config(menu=menu)
        return


    #
    # Menu actions
    #

    def new_project(self):
        return

    def open_project(self):
        return

    def quit(self):
        import os
        os._exit(0)
        return


if __name__=='__main__':
    Main()
