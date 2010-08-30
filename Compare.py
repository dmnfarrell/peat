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

"""Comparator class for comparing multiple proteins or mutants"""

import sys,os
from Tkinter import *
import Pmw
from GUI_helper import *
import math,re
from Ekin.Pylab import Options
from Ekin.Ekin_main import PlotPanel
from Actions import DBActions
from Yasara import YasaraControl

class Comparator(Frame, GUI_help):

    def __init__(self, parent=None, parentframe=None, DB=None, mode=None):
        """Instantiate class and check presence of DB data"""
        self.parent = parent
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
            self.main.title('PEAT Comparator')
            if DB == None:
                DB = self.openDB()
                if DB==None: return
        elif parentframe:
            self.main = Frame(master=parentframe)
            self.main.pack(fill=BOTH)
        else:
            self.main=Toplevel()
            self.main.title('Select Dataset')
            self.main.geometry('+300+200')
        if not self.parent:
            from Prefs import Preferences
            self.preferences=Preferences('PEAT',{})
        else:
            self.preferences=self.parent.preferences        
        self.DB = DB
        self.Opt = Options(redraw=None)
        self.doGUI()
        return

    def openDB(self):
        import tkFileDialog
        filename=tkFileDialog.askopenfilename(defaultextension='.fs',
                                       initialdir=os.getcwd(),
                                       filetypes=[("zodb fs","*.fs"),("All files","*.*")])
        if filename != None:        
            from Base import PDatabase
            DB = PDatabase(local=filename)
            return DB
        return None
        
    def doGUI(self):
        """Data selection dialog"""
        import tkFont
        dwin = self.main
        
        self.recs = self.DB.getRecs()
        fields = self.DB.getFields() + ['Structure']

        yscrollbar=Scrollbar(dwin,orient='vertical',width=12)
        yscrollbar.grid(row=1,column=2,sticky='news',padx=2)
        recsbox=Listbox(dwin,bg='white',
                         exportselection=0,
                         height=18,
                         width=20,
                         yscrollcommand=yscrollbar.set,
                         selectmode=EXTENDED)
        yscrollbar.config(command=recsbox.yview)
        Label(dwin,text='Records:').grid(row=0,column=0,sticky='news')
        recsbox.grid(row=1,column=0,columnspan=2,sticky='news',padx=1,pady=3)
        recsbox.config(state=NORMAL)

        y1scrollbar=Scrollbar(dwin,orient='vertical',width=12)
        y1scrollbar.grid(row=1,column=5,sticky='news',padx=2)
        colsbox=Listbox(dwin,bg='white',
                         exportselection=0,
                         height=20,
                         width=20,
                         yscrollcommand=y1scrollbar.set,
                         selectmode=EXTENDED)

        y1scrollbar.config(command=colsbox.yview)
        Label(dwin,text='Fields:').grid(row=0,column=3,sticky='news')
        colsbox.grid(row=1,column=3,columnspan=2,sticky='news',padx=1,pady=3)
        colsbox.config(state=NORMAL)
        dwin.rowconfigure(1, weight=1)
        dwin.columnconfigure(0, weight=1)
        for r in self.recs:
            recsbox.insert('end', r)
        for f in fields:
            colsbox.insert('end', f)

        def dorecslist():
            recsbox.delete(0, END)
            if srecsvar.get() == True:
                self.recs = self.parent.table.get_selectedRecordNames()
            else:
                self.recs = self.DB.getRecs()
            for r in self.recs:
                recsbox.insert('end', r)
            return

        other = Frame(dwin)
        other.grid(row=2,column=0,columnspan=7,sticky='news',padx=1,pady=3)
        srecsvar = BooleanVar(); srecsvar.set(False)
        Checkbutton(other, text='use only selected records', variable=srecsvar,
                        command=dorecslist).pack(side=LEFT,fill=BOTH)
        from Dialogs import PEATDialog
        self.pb=PEATDialog.createBusyBar(other)
        self.pb.pack(side=RIGHT,fill=BOTH,padx=5)
        d1 = Frame(dwin)
        d1.grid(row=1,column=6,columnspan=2,sticky='news',padx=1,pady=3)

        def display():
            srecs = [self.recs[int(item)] for item in recsbox.curselection()]
            sfields = [fields[int(item)] for item in colsbox.curselection()]
            slabel = entry1.getvalue()
            if len(srecs)==0 or len(sfields)==0:
                return
            g=groupby.getcurselection()

            if 'Structure' in sfields:
                self.displayStructure(srecs)         
                sfields.remove('Structure')    
            self.displayEkin(srecs, sfields, slabel, html=htmlvar.get(),
                               overlay=overlayvar.get(), groupby=g)
            
            return

        def preview():
            srecs = [self.recs[int(item)] for item in recsbox.curselection()]
            sfields = [fields[int(item)] for item in colsbox.curselection()]
            if len(srecs)==0 or len(sfields)==0:
                st.insert(END, 'Choose recs/fields!')
                return
            st.delete(1.0, END)
            slabel = entry1.getvalue()            
            sl = self.findNames(srecs, sfields, slabel)
            st.insert(END, str(sl))
            return

        def closedialog():
            dwin.destroy()
            if not self.parent:
                sys.exit()
            return

        entry1 = Pmw.EntryField(d1,
                    labelpos = 'n',
                    value = '',
                    command=preview,
                    label_text = 'Data label:')
        entry1.pack(side=TOP,fill=BOTH)
        Button(d1, text='Preview', command=preview).pack(side=TOP,fill=BOTH)
        Button(d1, text='Plot Options', command=self.Opt.plotSetup).pack(side=TOP,fill=BOTH)
        groupby = Pmw.OptionMenu (d1,
                labelpos = 'w',
                label_text = 'Group By:',
                items = ('field','record','none'),
                menubutton_width = 8 )
        groupby.pack(side=TOP,fill=BOTH)
        htmlvar = BooleanVar(); htmlvar.set(False)
        overlayvar = BooleanVar(); overlayvar.set(False)
        self.orientvar = BooleanVar(); self.orientvar.set(False)
        Button(d1, text='Display', command=display, bg='#FFFF99').pack(side=TOP,fill=BOTH)
        #Checkbutton(d1, text='send to web page', variable=htmlvar).pack(side=TOP,fill=BOTH)        
        Checkbutton(d1, text='overlay plots', variable=overlayvar).pack(side=TOP,fill=BOTH)
        Checkbutton(d1, text='orient vertical', variable=self.orientvar).pack(side=TOP,fill=BOTH)
        st = Pmw.ScrolledText(d1,
                labelpos = 'n',
                label_text='Data Labels found:',
                usehullsize = 1,
                hull_width = 200,
                hull_height = 220,
                text_wrap='char')
        st.pack(side=TOP,fill=BOTH)
        Button(d1, text='Close', command=closedialog).pack(side=TOP,fill=BOTH)

        return

    def findNames(self, srecs, sfields, slabel):
        """Get datasets list from search string"""
        if slabel == '':
            s = re.compile('[a-z]*\d',re.IGNORECASE)
        else:
            s = re.compile(slabel, re.IGNORECASE)
        slabels = []
        for r in srecs:
            rec = self.DB[r]
            for f in sfields:
                if f == 'Structure':
                    if rec.Structure != None:                    
                        slabels.append('has structure')
                    continue
                self.pb.update()
                ftype = self.DB.userfields[f]['field_type']

                if not rec.has_key(f):
                    continue
                E = rec[f]
                for d in E.datasets:
                    self.pb.update()
                    if s.match(d) and not d in slabels:
                        slabels.append(d)
        return slabels

    def displayStructure(self, srecs):
        """Display a structure"""
        app = self.preferences.get('molgraphApplication')
        path = self.preferences.get('molgraphAppPath')
        DBActions.DB = self.DB
        if DBActions.yasara != None:
            DBActions.yasara.Clear()

        colors = ['green',250,270,290,310,320,340]
        i=0
        for r in srecs:
            color = colors[i]
            i+=1
            DBActions.displayStructure(r, 'Structure', app, path, color, False)

        #try to align objects
        Y = DBActions.yasara
        if hasattr(Y, 'objects') and len(Y.objects)>1:
            Y.AlignMultiAll()
        return
   
    def displayEkin(self, srecs, sfields, slabel='', html=False,
                        overlay=False, cols=None, groupby='field'):
        """Format and display ekin plots from given recs/fields/label,
           plots can be grouped according to field if required"""

        from Ekin.Base import EkinProject
        tmppath = os.getcwd()
        slabels = self.findNames(srecs, sfields, slabel)
        all_labels = self.findNames(srecs, sfields, '')
        if len(slabels) == 0:
            return

        def insert(En, r, f, grp='field'):
            """reused in both loops"""
            E = self.DB[r][f]
            for d in all_labels:
                edata = E.getDataset(d)
                if edata != None:
                    if grp=='field':
                        name = r+' '+d
                    elif grp=='record':
                        name = f+' '+d
                    else:
                        name = r+' '+f+' '+d
                    En.insertDataset(edata, name,
                                     fit=E.getFitData(d), replace=True)

        #create ekinprojects grouped by field or record, or just one big one
        ekprojects = {}
        
        if groupby == 'none':
            En = ekprojects['all'] = EkinProject()
        if groupby == 'field' or groupby == 'none':
            for f in sfields:
                if groupby == 'field':
                    En = ekprojects[f] = EkinProject()
                for r in srecs:
                    if not self.DB[r].has_key(f):
                        continue
                    insert(En, r, f, grp=groupby)
        elif groupby == 'record':
            for r in srecs:
                En = ekprojects[r] = EkinProject()
                for f in sfields:
                    if not self.DB[r].has_key(f):
                        continue
                    insert(En, r, f, grp=groupby)

        if html == True:
            from Ekin.Web import EkinWeb
            ew = EkinWeb()
            dirs=['plots', 'csv']
            for p in dirs:
                if not os.path.exists(os.path.join(tmppath, p)):
                    os.mkdir(os.path.join(tmppath, p))
                    ew.showEkinPlots(project=En,
                              outfile=os.path.join(tmppath, 'plots.html'),
                              imgpath=os.path.join(tmppath, 'plots'),
                              path='plots/',
                              title='PEAT plots',
                              columns=3)
            import webbrowser
            webbrowser.open('plots.html',autoraise=1)
        else:
            if groupby=='field': t=sfields
            else: t=srecs
            if self.orientvar.get() == 1: orient=VERTICAL
            else: orient = HORIZONTAL
            cwin = ComparatorWin('Comparing %s' %t, labels=all_labels,
                                 selected=slabels, orient=orient)
            cwin.Opt = self.Opt
            for e in ekprojects:
                En = ekprojects[e]
                En.name = e                
                plotframe = cwin.add(E=En, selected=slabels, overlay=overlay)
            cwin.configSize()

        return       

class StructureControls(Frame):
    """A yasara controller, specifically for handling mutation views"""
    def __init__(self, parent, yasara):        
        Frame.__init__(self, parent, height=200, width=160)
        self.yasara = yasara
        self.showlabels=IntVar(); self.showlabels.set(1)               
        self.cbfr = Pmw.ScrolledFrame(self,                  
                    usehullsize = 1,
                    hull_width = 160,
                    hull_height = 100) 
        self.cbfr.grid(row=0,column=0,columnspan=2,sticky='news',padx=2,pady=2)   
        Label(self, text='Show Labels').grid(row=1,column=0,padx=2,pady=2)        
        c=Checkbutton(self, variable=self.showlabels, 
                         command=self.labelResidues)
        c.grid(row=1,column=1,sticky='news',padx=2,pady=2)
        c=Button(self,text='Align', command=self.align)
        c.grid(row=2,column=0,columnspan=2,sticky='news',padx=2,pady=2)                
        return

    def createResCheckbuttons(self, resnums):
        """Check buttons for residue list passed in ResiduesfromNames"""
        def updatelabels(tag, state):
            resnums=list(self.checkbuttons.getcurselection())   
            self.labelResidues(resnums)
            
        if hasattr(self, 'checkbuttons'):
            self.checkbuttons.destroy()
        self.checkbuttons = Pmw.RadioSelect(self.cbfr.interior(),
                buttontype = 'checkbutton',
                orient = 'vertical',
                labelpos = 'w',
                command = updatelabels,
                label_text = 'Residues:')       
        self.checkbuttons.pack(fill=BOTH,expand=1)
        for r in resnums:
            self.checkbuttons.add(r)            
        return
        
    def align(self):
        """try to align objects"""
        Y = self.yasara
        Y.AlignMultiAll()
        return
    
    def labelResidues(self, resnums=None):
        """Label residues from a list"""
        if resnums == None:
            resnums = self.resnums
        self.yasara.UnlabelRes('ALL')    
        if self.showlabels.get() == 0:
            #self.yasara.UnlabelRes('ALL')
            return
        for r in resnums:           
            try:                
                self.yasara.LabelRes(r,str(r),color='white')
            except:
                pass
        return

    def getResiduesfromNames(self, labels):
        """Given dataset names, extract residue numbers"""
        resnums=[]
        s = re.compile('\D',re.IGNORECASE)
        for l in labels:
            for n in s.split(l):
                if n != '' and n not in resnums:
                    resnums.append(n)  
        self.resnums = resnums 
        self.createResCheckbuttons(resnums)
        return resnums

    def markLabels(self, labels):
        """Use labels to mark residues"""
        resnums = self.getResiduesfromNames(labels)
        self.labelResidues(resnums)
        return
    
class DatasetSelector(Frame):
    """A dataset selection control"""
    def __init__(self, parent, labels, callback=None):        
        Frame.__init__(self, parent, width=160) 
               
        self.labels = labels
        self.labelvars = {}

        self.optmenu = Pmw.OptionMenu (self,
                labelpos = 'w',
                label_text = 'Datasets:',
                items = labels,
                menubutton_width = 10,
                command = callback)
        self.optmenu.setvalue('')
        self.optmenu.grid(row=0, column=0, columnspan=2, sticky='news')        
        
        datasetlist = Pmw.ScrolledFrame(self,
                    labelpos = 'n', label_text = 'Choose Datasets',
                    usehullsize = 1,
                    hull_width = 160,
                    hull_height = 350)        
        datasetlist.grid(row=1, column=0, columnspan=2, sticky='news')
        fr=datasetlist.interior()
        r=0; c=0; i=0
        for name in labels:
            if r>20:
                r=0
                c=c+2
            r=r+1
            Label(fr,text=name).grid(row=r,column=c)
            self.labelvars[name] = IntVar()
            self.labelvars[name].set(0)
            cb=Checkbutton(fr, onvalue=1,offvalue=0,variable=self.labelvars[name])
            cb.grid(row=r,column=c+1,padx=2,pady=2)                     
        b = Button(self, text="Apply", command=callback)
        b.grid(row=2,column=0,sticky='news',padx=2,pady=2)
        c=Button(self,text='Close', command=self.close)
        c.grid(row=2,column=1,sticky='news',padx=2,pady=2)
        c=Button(self,text='Select All', command=self.selectAll)
        c.grid(row=3,column=0,sticky='news',padx=2,pady=2)
        c=Button(self,text='Select None', command=self.selectNone)
        c.grid(row=3,column=1,sticky='news',padx=2,pady=2)
        self.columnconfigure(0,weight=1)
        return 

    def selectAll(self):
        for l in self.labelvars:
            self.labelvars[l].set(1)
        return
    
    def selectNone(self):
        for l in self.labelvars:
            self.labelvars[l].set(0)
        return
    
    def close(self):
        if self:
            self.destroy()
        return
        
    def getSelected(self):
        selected=[]
        for l in self.labels:
            if self.labelvars[l].get() == 1:
                selected.append(l)
        selected.append(self.optmenu.getvalue())        
        return selected     

    def setSelected(self, selected):        
        for l in selected:
            if l in self.labelvars:
                self.labelvars[l].set(1)                          
        return
    

class ComparatorWin(Frame):
    """Class to show multi plots"""
    def __init__(self, title, labels=None, selected=None,
                    orient=HORIZONTAL):
        Frame.__init__(self)
        self.main = Toplevel()
        self.main.title(title)
        self.width = 800
        self.main.geometry('%sx600+200+50' %self.width)
        self.orient = orient
        #add a selector for all the datasets in all plotframes, so we can replot subsets
        self.controls = Frame(self.main)
        #self.pw.add(self.controls)
        self.controls.pack(side=LEFT,fill=BOTH, expand=0)
        self.selector = DatasetSelector(self.controls, labels, callback=self.replot)
        self.selector.pack(side=TOP,fill=BOTH,expand=1)
        self.selector.setSelected(selected) 
        
        self.pw = PanedWindow(self.main,
                           orient=orient,
                           sashwidth=3,
                           showhandle=True,
                           opaqueresize=False)
        self.pw.pack(side=LEFT,fill=BOTH, expand=1)
        self.plotframes = {}
        if DBActions.yasara != None:
            self.yas = StructureControls(self.controls, DBActions.yasara)
            self.yas.markLabels(selected)
            self.yas.pack(side=TOP)
        else:
            self.yas=None
        return

    def add(self, E, selected=None, overlay=False):
        """Add ekin proj to plot"""
        
        total = len(E.datasets)
        if total == 0:
            return
        if selected == None:
            datasets = E.datasets
        else:
            datasets = self.getNames(selected, E)
            
        if len(datasets)>25:
            print 'too many plots'
            datasets = datasets[:25]
        elif len(E.datasets)==1:
            datasets = E.datasets[0]
        if len(datasets) == 0:
            datasets = E.datasets

        if overlay == True:
            plotopt = 3
            self.Opt.opts['legend'] = True
        else: plotopt = 2

        fr=Frame(self.pw)
        self.pw.add(fr)
        plotframe = PlotPanel(parent=fr, side=BOTTOM, tools=True)
        plotframe.setProject(E)
        plotframe.Opts.opts = self.Opt.opts            
        plotframe.Opts.opts['title'] = E.name
        plotframe.Opts.opts['normalise']=overlay
        plotframe.plotCurrent(datasets=datasets, plotoption=plotopt)
        self.plotframes[E.name] = plotframe       
        return plotframe

    def replot(self, event=None, selected=None):
        """Replot all frames with current datasets"""
        if selected == None:
            selected = self.selector.getSelected()         
        if len(selected) == 0:
            return
        for p in self.plotframes:
            pf = self.plotframes[p]
            E = pf.E
            #get all datasets matching values
            datasets = self.getNames(selected, E)            
            pf.Opts.opts['title'] = E.name
            pf.plotCurrent(datasets=datasets)
        if self.yas != None:    
            self.yas.markLabels(selected)
        return

    def getNames(self, selected, E):
        names=[]
        for d in E.datasets:
            for label in selected:
                if d.endswith(label) and label!='' and d not in names:
                    #print label, d
                    names.append(d)            
        return names
    
    def configSize(self):
        """Handle resize"""
        panes = self.pw.panes()
        #print len(panes)
        if self.orient == HORIZONTAL:
            size = float(self.main.winfo_width())/len(panes)
        else:
            size = float(self.main.winfo_height())/len(panes)

        i=0
        for fr in panes:          
            x=int(size*(i+1))
            self.pw.forget(fr)
            self.pw.paneconfig(fr, minsize=size)
            i+=1
        return


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="file",
                        help="Open a local db")
    opts, remainder = parser.parse_args()
    if opts.file != None and os.path.exists(opts.file):
        from PEATDB.Base import PDatabase
        DB = PDatabase(local=opts.file)
        C = Comparator(DB=DB)
        C.mainloop()
    else:
        C = Comparator()
        C.mainloop()
         
if __name__ == '__main__':
    main()

