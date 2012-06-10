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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

class SequenceLoader(Frame):
    """Dialog to load multiple sequences"""
    def __init__(self, parent, parentapp=None):
        self.parent = parent
        self.parentapp = parentapp
        Frame.__init__(self, parent)
        self.main=self
        self.makeGUI()
        #self.basesequence = self.parentapp.P.DNAseq
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
        self.m = PanedWindow(self.main,orient=HORIZONTAL,
                           sashwidth=2,
                           showhandle=True)
        self.m.pack(side=TOP,fill=BOTH,expand=1)
        self.filelist = Pmw.ScrolledListBox(self.m)
        self.filelist.component('listbox').configure(selectmode=EXTENDED)
        self.seqlist = Pmw.ScrolledListBox(self.m,
                            selectioncommand=self.update)
        self.seqlist.component('listbox').configure(selectmode=EXTENDED)        
        self.m.add(self.filelist)
        self.m.add(self.seqlist)
        btn=Button(self.main,text='Clear sequences',command=self.clear)
        btn.pack(side=TOP,fill=X)
        return

    def update(self, evt=None):
        #get sequences here and load them into canvas
        names = self.seqlist.getcurselection()
        seqs = [self.sequences[i] for i in names]
        sc = self.parentapp.sc
        sc.addSequences(seqs)
        sc.showMultiple()
        return

    def checkSequences(self, files):
        """Check if files hold sequences"""
        self.sequences = {}
        recs = []
        for f in files:
            handle = open(f, "rU")
            for record in SeqIO.parse(handle, "fasta"):
                #print record
                #self.sequences[record.id] = str(record.seq)
                recs.append(record)
            handle.close()
        baseseq = self.parentapp.P.DNAseq
        baserec = SeqRecord(Seq(baseseq), id='base')
        recs.append(baserec)
        print recs
        self.sequences = self.alignSequences(recs)
        
        return

    def alignSequences(self, recs):
        SeqIO.write(recs, 'temp.faa', "fasta")        
        align = clustalAlignment('temp.faa')
        print align
        seqs = {}
        for rec in align:
            seqs[rec.id] = str(rec.seq)
        return seqs

    def loadFiles(self):
        """Load multiple files"""
        filelist = Dialogs.openFilenames(self,ext=['txt','fasta','seq'])
        if type(filelist) == types.UnicodeType:
            #master_win is reference to parent root
            filelist = self.parentapp.master.tk.splitlist(filelist)
        self.checkSequences(filelist)
        self.seqlist.setlist(self.sequences.keys())
        self.filelist.setlist(filelist)
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
        self.seqlist.setlist(seqlist)       
        return

    def removeAll(self):
        self.seqlist.setlist([])
        return

    def clear(self):
        """Clear current displayed seq"""
        sc = self.parentapp.sc
        sc.clearMultiple()       
        return

def writeFastaFile(sequences):
    """Write sequences to fasta"""
    recs=[]
    for s in sequences:
        recs.append(SeqRecord(Seq(s)))
    SeqIO.write(recs, "test.faa", "fasta")
    return

def clustalAlignment(filename):
    name = os.path.splitext(filename)[0]
    from Bio.Align.Applications import ClustalwCommandline
    from Bio import AlignIO
    cline = ClustalwCommandline("clustalw", infile=filename)
    print 'performing clustal alignment..'    
    stdout, stderr = cline()
    align = AlignIO.read(name+'.aln', "clustal")
    return align

def PhyloTree():    
    from Bio import Phylo
    tree = Phylo.read("test.dnd", "newick")
    Phylo.draw_ascii(tree)

