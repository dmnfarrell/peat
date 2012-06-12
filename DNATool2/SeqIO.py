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
from Base import SequenceCanvas
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

class SequenceAlignmentTool(Frame):
    """App to align and handle multiple sequences"""
    def __init__(self, parent, parentapp=None):
        self.parent = parent
        self.parentapp = parentapp
        Frame.__init__(self, parent)
        self.main=self
        self.align = None
        self.makeGUI()
        #self.basesequence = self.parentapp.P.DNAseq
        return

    def menuBar(self):
        menuBar = Pmw.MenuBar(self.main,
                            hull_borderwidth = 1,
                            hull_width=300)
        menuBar.pack(side=TOP,fill=X)
        menuBar.addmenu('File', 'IO')
        menuBar.addmenuitem('File', 'command', 'Load files',
                        command = self.loadFiles,
                        label = 'Load files')
        menuBar.addmenuitem('File', 'command', 'Load zip',
                        command = self.loadZip,
                        label = 'Load zip')
        menuBar.addmenu('Alignment', '')
        menuBar.addmenuitem('Alignment', 'command', 'Show alignment',
                        command = self.showAlignmentResults,
                        label = 'Show last alignment results')
        return

    def makeGUI(self):
        self.menuBar()
        fr=Frame(self.main)
        fr.pack(side=TOP,fill=X)
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
        #print recs
        self.sequences = self.alignSequences(recs)
        return

    def getClustalPath(self):
        from Prefs import Preferences
        preferences = Preferences('DNATool2',{})        
        return preferences.get('clustalpath')

    def alignSequences(self, recs):
        SeqIO.write(recs, 'temp.faa', "fasta")
        command = self.getClustalPath()
        align = clustalAlignment('temp.faa', command)
        self.align = align
        seqs = {}
        for rec in align:
            seqs[rec.id] = str(rec.seq)
        return seqs

    def showAlignmentResults(self):
        print self.align.format("clustal")
        fr = Toplevel(self.main)
        at = AlignmentText(fr, labelpos = 'n',
                            label_text='Results',
                            text_wrap='word',
                            usehullsize = 1,
                            hull_width = 500,
                            hull_height = 300,)
        fr.geometry('600x600+200+200')
        at.pack(side=TOP,fill=BOTH,expand=1)
        at.settext(self.align.format("clustal"))
        at.formatText()
        return

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

class SequenceAnalysis(Frame):
    """Simple sequence pasing"""
    def __init__(self, parent, parentapp=None):
        self.parent = parent
        self.parentapp = parentapp
        Frame.__init__(self, parent)
        self.main=self
        self.getParentSequence()
        self.makeGUI()
        return

    def makeGUI(self):
        fr=Frame(self.main)
        fr.pack(side=TOP,fill=X)
        btn=Button(fr,text='Show Summary',command=self.showSummary)
        btn.pack(side=TOP,fill=X,expand=1,padx=1)
        return

    def getParentSequence(self):
        self.baseseq = self.parentapp.P.DNAseq
        return

    def showSummary(self):
        import pylab as plt
        from Bio.SeqUtils import GC
        print GC(self.baseseq)
        return

class AlignmentText(Pmw.ScrolledText):
    def __init__(self, parent=None, **kwargs):
        Pmw.ScrolledText.__init__(self, parent, **kwargs)
        return

    def formatText(self):
        self.tag_config('a', foreground="blue")
        self.tag_config('g', foreground="red")
        self.tag_config('t', background="yellow")
        self.tag_config('c', background="green")
        
        lines = int(self.index('end').split('.')[0]) - 1
        for l in range(lines):
            text = self.get(str(l)+'.0', str(l)+'.end')
            print text.split('|')
        return

class BlastInterface(Frame):
    """Simple interface to allow user to blast their sequence"""
    def __init__(self, parent, parentapp=None):
        self.parent = parent
        self.parentapp = parentapp
        Frame.__init__(self, parent)
        self.main=self
        self.makeGUI()
        self.getParentSequence()
        return

    def makeGUI(self):
        fr=Frame(self.main)
        fr.pack(side=TOP,fill=X)
        btn=Button(fr,text='Load Seq',command=self.loadSequence)
        btn.pack(side=LEFT,fill=X,expand=1,padx=1)
        btn=Button(fr,text='Save results',command=self.saveResults)
        btn.pack(side=LEFT,fill=X,expand=1,padx=1)
        btn=Button(fr,text='Run Blast',command=self.run)
        btn.pack(side=LEFT,fill=X,expand=1,padx=1)
        self.results = Pmw.ScrolledText(self.main,
                            labelpos = 'n',
                            label_text='Results',
                            text_wrap='word',
                            usehullsize = 1,
                            hull_width = 500,
                            hull_height = 300,)
        self.results.pack(side=TOP,fill=BOTH,expand=1,padx=1)
        btn=Button(self.main,text='Graphical summary',command=self.showGraphicalView)
        btn.pack(side=TOP,fill=X,padx=1)
        return

    def getParentSequence(self):
        self.baseseq = self.parentapp.P.DNAseq
        return

    def loadSequence(self):
        filename = Dialogs.openFilename(parent=self, ext='fasta')
        if not filename:
            return
        rec = readSequenceFile(filename)[0]
        self.baseseq = str(rec.seq)
        return

    def saveResults(self):

        return

    def loadBlastResults(self):

        return

    def run(self):

        baserec = SeqRecord(Seq(self.baseseq), id='base')
        from Bio.Blast import NCBIWWW

        self.pb = Dialogs.WaitDialog(self.main, message='Doing online BLAST. Please wait..',
                                 cancel=self.stopCurrent)
        self.pb.after(100, self.pb.updateValue())
        handle = NCBIWWW.qblast("blastn", "nt", baserec.seq)
        f = open("temp_blast.xml", "w")
        f.write(handle.read())
        f.close()
        handle.close()
        self.pb.close()

        #re-open handle from local file
        handle = open("temp_blast.xml")
        from Bio.Blast import NCBIXML
        blastrecs = NCBIXML.parse(handle)
        rec = blastrecs.next()
        self.alignedseqs = alignedseqs = []

        self.results.delete(1.0,END)
        self.insertAlignmentText(rec.alignments, self.results)
        self.alignments =rec.alignments
        return

    def insertAlignmentText(self, alignments, text):
        """Show formatted alignment info in text widget"""
        text.tag_config('a', foreground="blue")
        text.tag_config('r', foreground="red")
        text.tag_config('b', background="lightyellow")
        text.tag_config('c', background="lightgreen")
        n=100
        for align in alignments:
           for hsp in align.hsps:
                if hsp.expect < 0.01:
                    text.insert(END,'**Alignment** '+ 'score: ' +str(hsp.score)+ '\n',
                                ('a'))
                    text.insert(END, 'sequence: %s\n' %align.title,('r'))
                    for i in range(0, len(hsp.query), n):
                        text.insert(END, hsp.query[i:i+n],('b'))
                        text.insert(END,'\n')
                        text.insert(END, hsp.match[i:i+n])
                        text.insert(END,'\n')
                        text.insert(END, hsp.sbjct[i:i+n] ,('c'))
                        text.insert(END,'\n')
        return

    def stopCurrent(self):
        return

    def showGraphicalView(self):
        """Show canvas view of alignment"""
        print self.alignedseqs        
        AlignmentPlot(self.alignments)

def AlignmentPlot(self, alignments):
    import pylab as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cm
    vals = []
    for align in alignments:
       for hsp in align.hsps:
            #print align.length, hsp.score, hsp.expect
            vals.append((align.title[:30]+'..', hsp.score, hsp.expect))

    vals = sorted(vals, key=lambda x: x[2])
    vals = zip(*vals[:50])
    print vals
    f=plt.figure(figsize=(8,10))
    ax = plt.Axes(f, [.3,.05,.6,.9])
    f.add_axes(ax)
    names = vals[0]
    evals = vals[1]
    y = range(len(names))
    patches = ax.barh(y, evals, align='center')
    ax.set_yticks(y)
    ax.set_yticklabels(names)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(8)
    # set the colors of each bar based on the weights
    for val, patch in zip(evals, patches):
        # We need to normalize the "weight" value between 0-1
        # to generate an actual color...
        color = cm.jet(float(val) / max(evals))
        print color
        patch.set_facecolor(color)
    #cb = f.colorbar(scalarMap)
    f.subplots_adjust(hspace=0.5)
    plt.show()
    return

def chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

def getFormattedBLASTResults():
    return text

def readSequenceFile(filename, format='fasta'):
    recs = []
    for seqrecord in SeqIO.parse(filename, format):
        recs.append(seqrecord)
        print seqrecord.id
        print repr(seqrecord.seq)
        print len(seqrecord)
        return recs

def writeFastaFile(sequences):
    """Write sequences to fasta"""
    recs=[]
    for s in sequences:
        recs.append(SeqRecord(Seq(s)))
    SeqIO.write(recs, "test.faa", "fasta")
    return

def clustalAlignment(filename, command="clustalw"):
    name = os.path.splitext(filename)[0]
    from Bio.Align.Applications import ClustalwCommandline
    from Bio import AlignIO
    cline = ClustalwCommandline(command, infile=filename)
    print 'performing clustal alignment..'
    stdout, stderr = cline()
    align = AlignIO.read(name+'.aln', "clustal")
    return align

def PhyloTree():
    from Bio import Phylo
    tree = Phylo.read("test.dnd", "newick")
    Phylo.draw_ascii(tree)

def runBlast(record):
    from Bio.Blast import NCBIWWW
    result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
    return

