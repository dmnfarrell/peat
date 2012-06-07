#!/usr/bin/env python
#
# DNATool - A program for DNA sequence manipulation
# Copyright (C) 2012- Damien Farrell & Jens Erik Nielsen
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
# Email: farrell.damien_at_gmail.com

"""Methods for opening and wroting sequence formats"""

import os, types, string
from Tkinter import *
import Pmw
import Dialogs

class SequenceLoader(Frame):
    """Dialog to load multiple sequences"""
    def __init__(self, parent, parentapp=None):
        self.parent = parent
        self.parentapp = parentapp
        Frame.__init__(self, parent)
        self.main=self
        self.makeGUI()
        return

    def makeGUI(self):
        fr=Frame(self.main) 
        fr.pack(side=TOP,fill=X)     
        btn=Button(fr,text='Load files',command=self.loadFiles)
        btn.pack(side=LEFT,fill=X,expand=1,padx=1)
        btn=Button(fr,text='Load zip',command=self.loadZip)
        btn.pack(side=LEFT,fill=X,expand=1,padx=1)
        btn=Button(fr,text='Remove all',command=self.removeAll)
        btn.pack(side=LEFT,fill=X,expand=1,padx=1)
        self.details = Pmw.ScrolledListBox(self.main,                            
                            selectioncommand=self.update)
        self.details.component('listbox').configure(selectmode=EXTENDED)
        self.details.pack(side=TOP,fill=BOTH,expand=1)
        btn=Button(self.main,text='Clear sequences',command=self.clear)
        btn.pack(side=TOP,fill=X)
        return

    def update(self, evt=None):
        #get sequences here and load them into canvas
        name = self.details.getcurselection()[0]
        seq = self.sequences[name]
        sc = self.parentapp.sc
        sc.addSequences([seq])
        sc.showMultiple()
        return

    def checkSequences(self, filelist):
        """Check if files hold sequences"""
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        self.sequences = {}
        for f in filelist:
            handle = open(f, "rU")
            for record in SeqIO.parse(handle, "fasta"):
                print record.seq
                self.sequences[record.id] = str(record.seq)
            handle.close()            
        return

    def loadFiles(self):
        """Load multiple files"""
        filelist = Dialogs.openFilenames(self,ext=['txt','fasta','seq'])
        if type(filelist) == types.UnicodeType:
            #master_win is reference to parent root
            filelist = self.parentapp.master.tk.splitlist(filelist)
        self.checkSequences(filelist)
        self.details.setlist(self.sequences.keys())
        return

    def loadZip(self):
        zfilename = Dialogs.openFilename(self,ext=['zip'])
        if not zfilename: return
        zfile = zipfile.ZipFile( zfilename, "r" )
        zfile.printdir()
        seqlist=[]
        for info in zfile.infolist():
            fname = info.filename
            if fname.endswith(".clipped"):
                seqlist.append(fname)
                # decompress each file's data
                #data = zfile.read(fname)
                #print data
                x = os.path.split(fname)
        self.details.setlist(seqlist)       
        return

    def removeAll(self):
        self.details.setlist([])
        return

    def clear(self):
        """Clear current displayed seq"""

        return

