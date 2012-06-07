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


"""Base classes for DNATool"""

from Tkinter import *
import tkFont
import Pmw
import os, sys, types, pickle
import platform
import random
import Utilities
import Mutation
from PEATDB.Prefs import Preferences
import Dialogs

class Project(object):
    """Project class for DNATool data"""
    def __init__(self, data=None):
        if data==None:
            self.data = {}
        else:
            self.data = data
        self.data['DNAseq'] = ''
        return

    def load(self, filename):
        """load a db"""
        fd=open(filename,'rb')
        self.data = pickle.load(fd)
        fd.close()
        self.filename = filename
        return

    def save(self):
        fd=open(self.filename,'w')
        pickle.dump(fd, self.data)
        return

    def __getattr__(self, key):
        return self.data[key]

class Sequence(object):
    """Class to represent a sequence"""

    def __init__(self):
        """ """
        self.data = {}

        return

class SequenceCanvas(Pmw.ScrolledCanvas):
    """Canvas for drawing all graphical elements, sequences and restrictions sites"""
    def __init__(self, parent=None, parentapp=None,
                    width=400, height=400, bg='white', **kwargs):
        Pmw.ScrolledCanvas.__init__(self, parent, canvasmargin=2)
        self.parent = parent
        self.parentapp = parentapp
        self.platform = platform.system()
        self.width=width
        self.height=height
        self.canvas = self.component('canvas')
        self.canvas.configure(bg='#f6f9f6',height=height, width=width)
        self.createBindings()

        self.defaultprefs = {'seqfont':'Courier', 'seqfontsize':13, 'fontstyle':'bold',
                              'basescale':15, 'restrfont':'Arial',
                              'backgrcolor':'#f6f9f6','seqcoloroption':2}
        fonts = ['Courier','Arial','Fixed','Tahoma','Times','Verdana']
        #we use the following to create the prefs dialog
        #list of tuples of the form name,type,default,label
        self.defaultopts = [('seqfont','menu',fonts,'Sequence font'),
                            ('seqfontsize','scale',(6,30), 'Sequence font size'),
                            ('fontstyle','menu',['normal','bold','italic'],'Font style'),
                            ('basescale','scale',(8,50), 'Base scale'),
                            ('restrfont','menu',fonts,'Restr site font'),
                            ('backgrcolor','color','#f6f9f6','Background color')]
        self.loadPrefs()
        self.baserow = self.height/2
        self.seqstart = 80
        self.siterowstart = 70
        self.orfoffset = 15
        self.yscale = 1.1

        #self.sequence = Utilities.createDummySequence(100)
        self.alignedsequences = []
        self.update()
        return

    def createBindings(self):
        """Bindings"""
        c = self.canvas
        c.bind("<Button-1>",self.handleLeftClick)
        c.bind("<B1-Motion>", self.handleLeftMotion)
        c.bind("<Button-3>",self.handleRightClick)
        c.bind_all("<Control-KeyPress-c>", self.copySequence)
        #c.bind("<ButtonRelease-3>",self.handleRightRelease)

        '''c.bind("<ButtonRelease-1>",self.handle_left_release)

        c.bind("<Shift-Button-1>", self.handle_left_shift_click)
        c.bind("<Shift-ButtonRelease-1>", self.handle_left_shift_release)
        c.bind("<B3-Motion>", self.handle_right_motion)

        c.bind_all("<Control-KeyPress-x>",self.cut_DNA)
        c.bind_all("<Control-KeyPress-v>",self.paste_DNA)
        c.bind_all("<Delete>",self.delete_DNA) '''
        return

    def loadPrefs(self):
        """Load prefs from a preferences instance"""
        self.preferences = Preferences('DNATool2',{})
        temp = {}
        for key in self.defaultprefs:
            if not key in self.preferences.prefs:
                self.preferences.set(key, self.defaultprefs[key])
            temp[key] = self.preferences.get(key)
        self.applyPrefs(temp)
        return

    def savePrefs(self):
        """Save prefs to a preferences instance"""
        self.preferences = Preferences('DNATool2',{})
        for key in self.defaultprefs:
            self.preferences.set(key, self.__dict__[key])
        return

    def applyPrefs(self, prefs=None):
        """Apply prefs from a given dict"""
        if prefs == None: prefs=self.defaultprefs
        for key in prefs:
            self.__dict__[key] = prefs[key]
        return

    def handleLeftClick(self, evt):
        """Handle left click"""
        if hasattr(self, 'rightmenu'):
            self.rightmenu.destroy()
        c = self.canvas
        if len(c.gettags(CURRENT))==0:
            return
        self.markPosition(evt)
        x=evt.x; y=evt.y
        items = c.find_closest(evt.x, evt.y, halo=14)
        x1,y1,x2,y2 = c.bbox(CURRENT)
        self.selecteditems = c.find_enclosed(x1-1, y1-1, x2+1, y2+1)
        return

    def handleLeftMotion(self, evt):
        c = self.canvas
        items = self.selecteditems
        current = c.gettags(CURRENT)
        if 'sitelabel' in current:
            self.dragMove(evt, items)
        elif 'sequence' in current:
            self.selectSequence(evt)
        return

    def selectSequence(self, evt):
        """Select sequence"""
        c = self.canvas
        x = self.oldx
        x2 = c.canvasx(evt.x)
        height = self.basescale/1.1
        y = self.baserow+height/2
        c.delete(c.find_withtag('seqselection'))
        rect = c.create_rectangle(x,y, x2, y-height, width=0,
            tag=('seqselection'),fill='yellow')
        x1,y1,x2,y2 = c.bbox(rect)
        items = c.find_overlapping(x1+2, y1, x2-2, y2)
        seq = range(self.getSeqPositionFromCoords(x),self.getSeqPositionFromCoords(x2))
        self.colorSequence() #reset colors first
        for i in items:
            if 'seqtext' in c.gettags(i):
                c.itemconfig(i,fill='red')
                c.lift(i)
        self.selectedrange = seq
        return

    def highlightSequence(self, sequence):
        """Highlight section of sequence in different color"""
        c = self.canvas
        height = self.basescale/1.1
        y = self.baserow+height/2
        for s in sequence:
            x = self.getBasePosition(s)
            items =c.find_overlapping(x-2, y-2, x+2, y+2)
            for i in items:
                c.itemconfig(i,fill='red')
        return

    def getSeqPositionFromCoords(self,x):
        """From X and Y coordinate, return the DNA base number"""
        #if abs(y-self.baserow)>5:
        #    return None
        pos = int(float(x-self.seqstart+4.0)/self.basescale)
        return pos

    def markPosition(self, evt):
        c = self.canvas
        self.oldx = c.canvasx(evt.x)
        self.oldy = c.canvasy(evt.y)
        self.selectedtags = c.gettags(CURRENT)
        return

    def dragMove(self, evt, items=None, tags=None):
        """Handle mouse drag"""
        c = self.canvas
        x = c.canvasx(evt.x)
        y = c.canvasy(evt.y)
        diffx = x - self.oldx
        diffy = y - self.oldy
        self.oldx = x
        self.oldy = y
        name = c.gettags(CURRENT)[0]
        for i in items:
            c.move(i, diffx, diffy)
            self.drawSiteConnector(x, y, name, move=True)
        return

    def handleRightClick(self, evt):
        self.rightmenu = self.popupMenu(evt)
        return

    def popupMenu(self, evt):
        """Right click shows popup menu"""
        popupmenu = Menu(self.canvas, tearoff = 0)
        def popupFocusOut(evt):
            popupmenu.unpost()
        popupmenu.add_command(label="Show Prefs", command= self.showPrefs)
        popupmenu.bind("<FocusOut>", popupFocusOut)
        popupmenu.focus_set()
        popupmenu.post(evt.x_root, evt.y_root)
        return popupmenu

    def addSequences(self, sequences):
        self.alignedsequences = sequences
        return

    def update(self, sequence=None, restrictionsites=None, primers=None):
        """Update display of the current sequence(s)"""

        #if sequence == None:
        #    return

        if self.parentapp != None:
            P = self.parentapp.P
            if P.DNAseq == '':
                return
        sequence = P.DNAseq
        enzymes = P.used_enzymes
        cutpos = P.cut_pos

        #remove previous
        self.clear()
        self.calculateScaling()
        self.drawDefaultLabels()
        #update sequence
        self.showSequence(sequence)
        #if self.ORFselected == True:
        self.showAA3Seq(sequence)
        #draw restr sites
        self.showRestrictionSites(enzymes, cutpos)
        #update primers if selected
        self.showPrimers()
        self.canvas.configure(bg=self.backgrcolor)
        return

    def clear(self):
        """Clear all"""
        c = self.canvas
        c.delete(ALL)
        return

    def calculateScaling(self):
        """Calculate font size dependant y-scaling to prevent crowding"""
        self.yscale = self.seqfontsize/10.0
        self.yscale=pow(self.yscale,0.05)
        return

    def getCurrentFont(self):
        """Generate current font from settings"""
        fontsize = self.seqfontsize
        fontstyle = self.fontstyle
        #create a font object
        if fontstyle != 'italic':
            seqfont = tkFont.Font (family=self.seqfont, size=fontsize,
                                    weight=fontstyle)
        else:
            seqfont = tkFont.Font (family=self.seqfont, size=fontsize,
                                    slant="italic")
        return seqfont

    def showSequence(self, sequence=None, row=None):
        """Show sequence"""
        if sequence == None:
            return
        if row == None:
            row = self.baserow
        c = self.canvas
        seqfont = self.getCurrentFont()

        #dna sequence
        seq = sequence
        print sequence
        for i in range(len(seq)):
            pos = self.getBasePosition(i)
            item = c.create_text(pos,row,
                            font=seqfont,text=seq[i],
                            fill='black',anchor='w',tag=('seqtext','sequence'))
        self.seqend = self.getBasePosition(len(seq))-self.seqstart
        self.colorSequence()
        self.drawSeqBackground()
        self.canvas.configure(scrollregion=(0,0,self.seqend+100,self.height))
        self.centerPage()
        return

    def showMultiple(self):
        """Show multiple sequences aligned to the base one"""
        print self.alignedsequences
        if len(self.alignedsequences) == 0:
            return
        row = self.baserow - 30
        for s in self.alignedsequences:
            self.showSequence(s, row)
        return
            
    def showAASeq(self):

        return

    def showAA3Seq(self, sequence=None):
        """Show 3 letter code amino acid sequence"""
        if sequence == None:
            return
        c = self.canvas
        seqfont = self.getCurrentFont()

        y = (self.baserow+self.orfoffset)*self.yscale
        frame = 0
        AAseqs3,AAseqs1 = Mutation.translate(sequence)
        aaseq = AAseqs3[frame]
        #print aaseq
        for i in range(len(aaseq)):
            pos = self.getAAPosition(i)#,frame+1)
            c.create_text(pos,y,
                            font=seqfont,text=aaseq[i],
                            fill='black',tag=('aasequence'))
        for i in range(len(aaseq)):
            if i%5!=0:
                continue
            pos = self.getAAPosition(i)
            c.create_text(pos,y+45,
                            font='Courier 11',text=i,
                            fill='black',tag=('aasequence'))
            c.create_line(pos,y+35,pos,y+20,fill='red',width=2)
        return

    def colorSequence(self):
        """Color the dna sequence"""
        c = self.canvas
        items = c.find_withtag('seqtext')
        option = 2
        count=0
        clr = '#66CC00'
        for i in items:
            if option == 1:
                c.itemconfig(i, fill=clr)
            elif option == 2:
                if count % 3 == 0 and clr == '#66CC00':
                    clr='#999900'
                elif count % 3 == 0 and clr=='#999900':
                    clr='#66CC00'
                c.itemconfig(i, fill=clr)
            count=count+1
        return

    def drawDefaultLabels(self):
        """Draw objects that are permanent"""
        c = self.canvas
        item = c.create_text(5,self.baserow,font='Courier 10',text='DNA seq',
                        fill='darkgreen',anchor='w')
        y = (self.baserow+self.orfoffset)*self.yscale
        item = c.create_text(5,y,font='Courier 10',text='ORF',
                        fill='black',anchor='w')
        item = c.create_text(5,y+45,font='Courier 10',text='AA no.',
                        fill='gray',anchor='w')
        return

    def drawSeqBackground(self, color='white'):
        """Draw a background for a sequence"""
        c = self.canvas
        x = self.seqstart
        height = self.basescale
        y = self.baserow+height/2
        end = self.seqend
        rect = c.create_rectangle(x,y, x+end, y-height,
                                   width=0,tag=('seqbackground','sequence'),fill=color)
        c.tag_lower(rect)
        return rect

    def drawRestrictionSite(self, name, positions):
        """Draw restriction site(s) for a single enzyme"""
        c = self.canvas
        font = restrfont = str(self.restrfont)+' '+str(self.seqfontsize-2)
        done=[]
        if len(positions) == 1:
            #font = font+' underline'
            fill = '#FF9999'
        else:
            fill = '#F6FCE1'
        for pos in positions:
            uid = name+str(pos[0])
            self.restrictionsites[uid] = pos[0]
            if pos[0] in done: continue
            done.append(pos[0])
            x = self.getBasePosition(pos[0])
            y = self.baserow-self.siterowstart
            text = self.create_text(x,y,text=name,tags=(uid,'sitelabel'),
                                        font=font,anchor='nw')

            x1,y1,x2,y2 = box = c.bbox(text)
            rect = c.create_rectangle(box,fill=fill,outline='gray',
                                        tags=(uid,'rect'))
            items = c.find_withtag(uid)
            overlapping = c.find_overlapping(x1,y1,x2,y2)[:-2]
            #print pos, items, overlapping
            inc=-10
            while len(overlapping) > 1:
                for i in items:
                    c.move(i,0,inc)
                x1,y1,x2,y2 = c.bbox(i)
                overlapping = c.find_overlapping(x1,y1,x2,y2)[:-2]
                y=y+inc
            c.lift(text)
            self.drawSiteConnector(x,y,uid)
        return text

    def drawSiteConnector(self,x2,y2,uid,move=False):
        """Draw line connecting site label to sequence"""
        c = self.canvas
        pos = self.restrictionsites[uid]
        x1 = self.getBasePosition(pos)
        y1 = self.baserow-15

        if move==True:
            old = list(set(c.find_withtag(uid))&set(c.find_withtag('line')))
            for i in old:
                c.delete(i)
        line = c.create_line(x1,y1,x1,y2, fill='blue',
                              width=2,stipple='gray25',
                              tag=(uid,'line'))
        c.tag_lower(line)
        line = c.create_line(x1,y2,x2,y2, fill='blue',
                              width=2,stipple='gray25',
                              tag=(uid,'line'))
        c.tag_lower(line)
        return

    def moveItemByName(self, name, y):
        c = self.canvas
        items = c.find_withtag(name)
        for i in items:
            c.move(i, 0, y)

        return

    def showRestrictionSites(self, enzymes=None, cutpos=None):
        """Draw all restriction sites"""
        if enzymes==None:
            return
        self.restrictionsites = {}
        c = self.canvas
        for e in enzymes:
            if cutpos.has_key(e):
                positions = cutpos[e]
                r = self.drawRestrictionSite(e, positions)
        print self.restrictionsites
        return

    def getBasePosition(self, pos):
        """Return the x position of the nth base on the canvas"""
        return self.seqstart + float(pos) * self.basescale

    def getAAPosition(self, pos, frame=None):
        pos = pos*3+1 #+float(frame)
        return self.getBasePosition(pos)

    def showPrimers(self):
        """Draw primers"""
        return

    def centerPage(self):
        """Center canvas on sequence after redraw"""
        top, bottom = self.yview()
        size = bottom - top
        if size==0: size=0.4
        middle = size * 0.4
        self.yview('moveto', middle)
        #print top, bottom
        #print middle
        return

    def zoomIn(self):
        """Zoom in by enlarging all elements"""
        c = self.canvas
        x1,x2 = self.xview()
        y1,y2 = self.yview()
        self.basescale = self.basescale*1.1
        self.seqfontsize = self.seqfontsize+1
        self.update()
        return

    def zoomOut(self):
        """Zoom out by reducing all elements"""
        if self.seqfontsize<7:
            return
        self.basescale = self.basescale/1.1
        self.seqfontsize = self.seqfontsize-1
        self.update()
        return

    def copySequence(self, evt=None):
        """Copy sequence text"""
        seqdata = ''.join([self.sequence[i] for i in self.selectedrange])
        self.parent.clipboard_clear()
        self.parent.clipboard_append(seqdata)
        return

    def showPrefs(self, evt=None):
        """Show preferences window"""
        fr = Toplevel(self.parent)
        fr.title('Sequence window preferences')
        def func():
            self.applyPrefs(self.prefsdialog.getValues())
            self.update()
        curr = {}
        for key in self.defaultprefs:
            curr[key] = self.__dict__[key]
        dlg = self.prefsdialog = PrefsDialog(fr, self, options=self.defaultopts,
                                    defaults=curr, callback=func)
        dlg.pack(fill=BOTH,expand=1)
        fr.geometry(dlg.geometry)
        fr.grab_set()
        fr.transient(self.parent)
        return

class ORFOverview(Frame):
    """CLass showing overview of ORF"""
    def __init__(self, parent, parentapp=None, sequence=None):
        self.parent = parent
        self.parentapp = parentapp
        Frame.__init__(self, parent)
        self.main = self
        self.font = 'Arial 9'
        self.xstart=80
        self.width=300
        #reference to current project
        if parentapp != None:
            self.P = self.parentapp.P
        self.makeGUI()
        self.bind('<Configure>', self.update)
        return

    def makeGUI(self):
        self.sc = Pmw.ScrolledCanvas(self.main)
        self.sc.pack(side=TOP,fill=BOTH,expand=1)
        self.canvas = self.sc.component('canvas')
        self.canvas.configure(bg='#f6f9f6',relief='groove',height=200,width=self.width)
        return

    def update(self, evt=None):
        """Draw all elements"""
        sequence = self.P.DNAseq
        #print sequence
        self.width=self.winfo_width()
        start = self.xstart
        end = self.width-20
        self.scale = (end-start)/float(len(sequence))
        c = self.canvas
        c.delete(ALL)
        import Mutation
        AAseqs3,AAseqs1 = Mutation.translate(sequence)
        #print AAseqs1
        i=1

        for seq in AAseqs1:
            row = 40+i*30
            c.create_text(start-60,row,
                          font=self.font,text='Frame '+str(i),
                          fill='black',anchor='w',tag=('label'))

            c.create_line(start,row,end,row,width=2,
                          fill='darkgreen',tag=('line'))
            pos=1
            for letter in seq:
                if letter=='*':
                    x = self.getAAPosition(pos)
                    print pos,x, self.scale, len(sequence)
                    c.create_line(x,row,x,row-10,fill='red',width=2)
                    c.create_text(x,row-20,text='%d' %(pos),
                                    font='Arial 7',anchor='n')
                pos+=1
            i+=1
        return

    def getAAPosition(self, pos, frame=None):
        pos = pos*3+1 #+float(frame)
        return self.getBasePosition(pos)

    def getBasePosition(self, pos):
        return self.xstart + float(pos) * self.scale

class PrefsDialog(Dialogs.GenericDialog):
    def __init__(self, parent, parentapp=None, options=[], defaults=None, callback=None):
        self.geometry = '300x300+500+250'
        self.parentapp = parentapp
        Dialogs.GenericDialog.__init__(self, parent, options, defaults, callback)
        b=Button(self.buttonframe, text='Save', command=self.parentapp.savePrefs)
        b.pack(side=LEFT,fill=X,expand=1,padx=1,pady=1)
        return

