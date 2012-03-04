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

import os
from Tkinter import *
import tkFileDialog
import Pmw
import fileinput

class TextEditor(Frame):
    """Class to provide a frame for text editing"""

    def __init__(self, parent=None, filename=None, title=None):
        self.parent=parent
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
        else:
            self.main=Toplevel()
            self.master=self.main
        self.main.title('Text Edit')
        ws = self.main.winfo_screenwidth()
        hs = self.main.winfo_screenheight()
        w=600; h=600
        x = (ws/2)-(w/2)+40; y = (hs/2)-(h/2)
        self.main.geometry('%dx%d+%d+%d' % (w,h,x,y))
        self.main.protocol('WM_DELETE_WINDOW',self.quit)
        self.makeGUI()
        return

    def makeGUI(self):
        self.text = Pmw.ScrolledText(self.main,
                labelpos = 'n',
                label_text='Text',
                usehullsize = 1,
                hull_width = 300,
                hull_height = 400,
                text_wrap='word')
        self.text.grid(row=0,column=0,columnspan=2,sticky='NEWS',padx=2,pady=2)
        self.main.columnconfigure(0,weight=1)
        self.main.rowconfigure(0,weight=1)
        buttonpanel = Frame(self.main)
        buttonpanel.grid(row=1,column=0,columnspan=2,sticky='NEWS',padx=2,pady=2)
        self.doButtons(buttonpanel)
        return

    def doButtons(self, frame):
        """Buttons"""
        buttons=[['Open',self.openFile],
                 ['Save',self.save],
                 ['Save As',self.saveas],
                 ['Close',self.quit]]

        for button,command in buttons:
            x=Button(frame,text=button,command=command)
            x.pack(side=LEFT,fill=BOTH,padx=2,pady=2)
        return

    def openFile(self, filename=None):
        """Displays text from a simple text file"""
        if filename == None:
            filename = self.openFilename()
        if not filename:
            return
        for line in fileinput.input(filename):
            self.text.insert(END, line)
        self.filename = filename
        self.main.title(filename)
        return

    def save(self, filename=None):
        if filename == None:
            filename = self.filename
        contents = self.text.get(1.0, END)
        fd=open(filename,'w')
        fd.write(contents)
        fd.close()
        self.filename = filename
        return

    def saveas(self):
        filename = self.saveFilename()
        if not filename:
            return
        self.save(filename)
        return

    def saveFilename(self, ext='.txt'):
        filename=tkFileDialog.asksaveasfilename(defaultextension=ext,
                                              initialdir=os.getcwd(),
                                              filetypes=[("txt files","*.txt"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        return filename

    def openFilename(self, ext='.txt'):
        filename=tkFileDialog.askopenfilename(defaultextension=ext,
                                              initialdir=os.getcwd(),
                                              filetypes=[("text files","*.txt"),
                                                         ("csv files","*.csv"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        return filename

    def quit(self,event=None):
        """Close the frame"""
        self.main.destroy()
        return

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="file",
                        help="Raw file", metavar="FILE")
    opts, remainder = parser.parse_args()
    app = TextEditor()

    if opts.file != None:
        app.openFile(opts.file)
    app.mainloop()

if __name__ == '__main__':
    main()
