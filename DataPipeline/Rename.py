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

'''Module for batch file renamine'''

import os, string
import re, glob
from Tkinter import *
import Pmw
import tkFileDialog

def doRename(files=None, wildcard=None, pattern='', replacement='', rename=False):
    newfiles = []
    if files==None:
        files = glob.glob(wildcard)
    for pathname in files:
        basename= os.path.basename(pathname)
        new_filename= re.sub(pattern, replacement, basename)
        if new_filename != basename:
            newfiles.append(new_filename)
            if rename == True:
                os.rename(pathname,
                          os.path.join(os.path.dirname(pathname),
                          new_filename))
    return newfiles

def doFindReplace(files=None, wildcard=None, find='', replace='', rename=False):
    newfiles = []
    if files==None:
        files = glob.glob(wildcard)
    for pathname in files:
        basename= os.path.basename(pathname)
        new_filename = basename.replace(find,replace)
        newfiles.append(new_filename)
        if new_filename != basename:

            if rename == True:
                os.rename(pathname,
                          os.path.join(os.path.dirname(pathname),
                          new_filename))
    return newfiles

def constructRegex(inputstr):
    return

class BatchRenameApp(Frame):
    """Batch renaming app"""
    def __init__(self, parent=None):
        self.parent=parent
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
        else:
            self.main=Toplevel()
            self.master=self.main
        self.main.title('Batch Rename')
        ws = self.main.winfo_screenwidth()
        hs = self.main.winfo_screenheight()
        w = 800; h=500
        x = (ws/2)-(w/2); y = (hs/2)-(h/2)
        self.main.geometry('%dx%d+%d+%d' % (w,h,x,y))
        self.doGUI()
        self.currentdir = os.path.expanduser('~')
        return

    def doGUI(self):
        self.m = PanedWindow(self.main,
                           orient=HORIZONTAL,
                           sashwidth=3,
                           showhandle=True)
        self.m.pack(side=TOP,fill=BOTH,expand=1)
        self.listbox = Pmw.ScrolledListBox(self.m,
                labelpos='n',
                label_text='File names',
                listbox_height = 8,
                usehullsize = 1,
                hull_width = 400)
        self.listbox.component('listbox').configure(selectmode=EXTENDED)
        self.m.add(self.listbox)
        self.preview = Pmw.ScrolledListBox(self.m,
                labelpos='n',
                label_text='New names',
                usehullsize = 1,
                hull_width = 200)
        self.m.add(self.preview)

        fr=Frame(self.main,pady=4)
        b=Button(fr,text='Add Folder',command=self.addFolder)
        b.pack(side=LEFT,fill=BOTH,)
        b=Button(fr,text='Clear',command=self.clear)
        b.pack(side=LEFT,fill=BOTH)
        self.filepattern = Pmw.EntryField(fr,
                labelpos = 'w',
                value = '*.*',
                label_text = 'Wildcard:',
                command = self.refresh)
        self.filepattern.pack(side=LEFT,fill=BOTH)
        self.filepattern.component("entry").configure(width=8)
        self.findtext = Pmw.EntryField(fr,
                labelpos = 'w',
                value = ' ',
                label_text = 'Find:',
                command = self.refresh)
        self.findtext.pack(side=LEFT,fill=BOTH)
        self.findtext.component("entry").configure(width=8)
        self.replacetext = Pmw.EntryField(fr,
                labelpos = 'w',
                value = '.',
                label_text = 'Replace With:',
                command = self.refresh)
        self.replacetext.pack(side=LEFT,fill=BOTH)
        self.replacetext.component("entry").configure(width=8)

        b=Button(fr,text='Preview',command=self.dopreview)
        b.pack(side=LEFT,fill=BOTH)
        b=Button(fr,text='Execute',bg='lightblue',command=self.execute)
        b.pack(side=LEFT,fill=BOTH)
        fr.pack(side=TOP,fill=BOTH)
        return

    def addFolder(self,path=None):
        if path==None:
            path = tkFileDialog.askdirectory(parent=self.main,
                                            initialdir=self.currentdir,
                                            title='Select folder')
        if path:
            self.path = path
            self.refresh()
            self.currentdir = path
        return

    def refresh(self):
        fp = self.filepattern.getvalue()
        flist = glob.glob(os.path.join(self.path,fp))
        self.listbox.setlist(flist)
        self.dopreview()
        return

    def dopreview(self):
        self.preview.delete(0,END)
        flist = self.listbox.get()
        find=self.findtext.getvalue()
        repl=self.replacetext.getvalue()
        new = doFindReplace(files=flist, find=find, replace=repl)
        for f in zip(flist,new):
            if f[0] != f[1]:
                self.preview.insert(END,f[1])
        return

    def clear(self):
        self.listbox.delete(0,END)
        self.preview.delete(0,END)
        self.path = None
        return

    def execute(self):
        #askyesno
        flist = self.listbox.get()
        find=self.findtext.getvalue()
        repl=self.replacetext.getvalue()
        doFindReplace(files=flist, find=find, replace=repl, rename=True)
        self.refresh()
        return

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-d", "--dir", dest="directory",
                        help="Folder of raw files")

    opts, remainder = parser.parse_args()
    app = BatchRenameApp()
    if opts.directory != None:
        app.addFolder(opts.directory)
    app.mainloop()

if __name__ == '__main__':
    main()
