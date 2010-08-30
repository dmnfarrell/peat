#!/bin/env python

#
# $Id: DataBase_GUI.py 18 2005-04-17 23:28:22Z nielsen $
#
# Database GUI - interface to the Database routines
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
