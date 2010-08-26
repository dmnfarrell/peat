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

class Main_pKa_Calc_Window(Frame):
    #
    # Class setting up the main window
    #
    
    def __init__(self):
        #
        # Define the look and feel
        #
        self.font="Times 12 bold"
        self.bg_colour="grey"
        self.fg_colour="black"
        #
        # do the window
        #
        self.data={'pdbfilename':None}
        self.do_window()
     
        return

    def do_window(self):
        
        #
        # Open the main window
        #
        Frame.__init__(self)
        self.master.title("pKa Design - Main Window")

        #self.master.geometry("400x150")
        # Set up the main frame
        #self.master.rowconfigure(0,weight=1)
        #self.master.columnconfigure(0,weight=1)
        #self.master.grid(sticky=W)

        #
        # Text box til top
        #
        label1=Label(self.master, text="PDB file",font=self.font)
        label1.grid(row=0,column=0, sticky=W)

        
            
        self.pdbfield=Entry(self.master)
        #
        # If we have a filename then insert it
        #
        if self.data['pdbfilename']:
            self.pdbfield.insert(INSERT,self.data['pdbfilename'])
        self.pdbfield.grid(row=0,column=1,sticky=W)

        #
        # Button for selecting the pdb file
        #
        self.loadbutton=Button(self.master,text='Browse',
                               command=self.pdb_browse,
                               font=self.font,
                               fg=self.fg_colour,
                               bg=self.bg_colour)
        self.loadbutton.grid(row=0,column=2,sticky=W)

        #
        # Button for loading the PDB file
        #
        loadbutton=Button(self.master,text='Load',command=self.loadpdb)
        loadbutton.grid(row=1,column=1,sticky=W)

        #
        # Exit button
        #
        exitbutton=Button(self.master,text='Exit',command=self.exit)
        exitbutton.grid(row=1,column=0,sticky=W)

        return
    

    def pdb_browse(self):
        #
        # Open the dialog for selecting the PDB file
        #
        pdbfilename=tkFileDialog.askopenfilename(defaultextension='.pdb',
                                                 initialdir=os.getcwd(),
                                                 filetypes=[("PDB file","*.pdb"),
                                                            ("All files","*.*")],
                                                 parent=self.master)
        return

    def loadpdb(self):
        self.master.destroy()
        X=pKa_setup(self.data)
        return

    def update(self):
        #
        # Redraw the main window
        #
        print 'Updating main window'
        #
        self.do_window()
        return

    def exit(self):
        #
        # Exit the program
        #
        os._exit(0)
        return

    


class pKa_setup(Frame):

    def __init__(self,data):
        self.do_window()
        return

    def do_window(self):
        Frame.__init__(self)
        self.master.title("pKa calculations: Setup")
        return
