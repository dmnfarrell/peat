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


from PEATDB.Ekin.Base import *
from PEATDB.Tables import TableCanvas
from PEATDB.TableModels import TableModel
from PEATDB.Ekin.Ekin_main import EkinApp, PlotPanel
from PEATDB.Ekin.Pylab import Options 
import os, sys, math, random, glob, numpy, string, types
from datetime import datetime
import ConfigParser, csv
from Tkinter import *
import Pmw
import tkFileDialog
from PEATDB.GUI_helper import *

class PipeApp(Frame, GUI_help):
    """Data pipe GUI for importing and fitting of raw data.
       This class uses ekin provided an automated pipeline for fitting
       raw text data files and propagating errors.
       Uses a config file to store the pipeline settings"""
                                  
    def __init__(self, parent=None, rawfile=None, conffile=None):
        self.parent=parent    
        self.tableformat = {'cellwidth':50, 'thefont':"Arial 10",
                            'rowheight':16, 'editable':False,
                            'rowselectedcolor':'yellow'}           
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
        else:
            self.main=Toplevel()
            self.master=self.main
        self.main.title('DataPipeline Desktop')
        self.main.geometry('800x600+200+100')
        self.main.protocol('WM_DELETE_WINDOW',self.quit)
        
        #pipeline object is used for everything except gui stuff
        self.p = Pipeline(conffile)
        self.setupVars()
        self.setupGUI()
        #redirect stdout to log win
        self.log.delete(1.0,END)
        sys.stdout = self
        #sys.stderr = self   
        return

    def setupVars(self):
        """tk vars"""
        self.conffilevar = StringVar()
        self.queuefilesvar = IntVar()
        self.currfilevar = StringVar()
        return
    
    def setupGUI(self):
        """Do GUI elements"""
        self.createMenuBar()
        self.m = PanedWindow(self.main,
                           orient=HORIZONTAL,
                           sashwidth=3,
                           showhandle=True)
        self.m1 = PanedWindow(self.m,
                           orient=VERTICAL,
                           sashwidth=3,
                           showhandle=True)
        self.m2 = PanedWindow(self.m,
                           orient=VERTICAL,
                           sashwidth=3,
                           showhandle=True)
        self.m.pack(side=TOP,fill=BOTH,expand=1) 
        self.m.add(self.m1)
        self.rawcontents = Pmw.ScrolledText(self.m1,
                labelpos = 'n',
                label_text='Raw File Contents',
                rowheader=1,
                columnheader=1,
                Header_foreground = 'blue',
                rowheader_width = 3,
                usehullsize = 1,
                hull_width = 500,
                hull_height = 300,
                text_wrap='none')
        self.m1.add(self.rawcontents)
        self.previewer = PlotPreviewer(app=self)
        self.m1.add(self.previewer)        
        
        self.m.add(self.m2)
        self.queueFrame = queueManager(app=self)
        self.m2.add(self.queueFrame)
        self.log = Pmw.ScrolledText(self.m2,
                labelpos = 'n',
                label_text='Logs',
                usehullsize = 1,
                hull_width = 400,
                hull_height = 500,
                text_wrap='word')       
        self.m2.add(self.log)
        self.infopane = Frame(self.main,height=20)
        self.infopane.pack(side=BOTTOM,fill=BOTH,pady=4)
        self.updateinfoPane()
        Label(self.infopane,text='Conf file:').pack(side=LEFT)
        Label(self.infopane,textvariable=self.conffilevar,fg='darkblue').pack(side=LEFT,padx=4)
        Label(self.infopane,text='Files in queue:').pack(side=LEFT,padx=4)
        Label(self.infopane,textvariable=self.queuefilesvar,fg='darkblue').pack(side=LEFT)
        Label(self.infopane,text='Current file:').pack(side=LEFT,padx=4)
        Label(self.infopane,textvariable=self.currfilevar,fg='darkblue').pack(side=LEFT)        
        return

    def updateinfoPane(self):
        self.conffilevar.set(self.p.conffile)
        self.queuefilesvar.set(len(self.p.queue))
        self.currfilevar.set(self.p.filename)
        return
    
    def createMenuBar(self):
        """Create the menu bar for the application. """
        self.menu=Menu(self.main)
        self.file_menu={'01Open Raw File(s)':{'cmd':self.openRaw},
                         '02Load Config File':{'cmd':self.loadConfig},
                         '03Edit Current Config':{'cmd':self.editConfig},
                         '04Create Config':{'cmd':self.createConfig},                        
                         '05Quit':{'cmd':self.quit}}                         
        self.file_menu=self.create_pulldown(self.menu,self.file_menu)
        self.menu.add_cascade(label='File',menu=self.file_menu['var'])
        self.presetsMenu()                
        self.run_menu={'01Execute':{'cmd': self.execute},                       
                        '02Open results in Ekin':{'cmd':self.openEkin}}       
        self.run_menu=self.create_pulldown(self.menu,self.run_menu)        
        self.menu.add_cascade(label='Run',menu=self.run_menu['var'])
        self.queue_menu={'01Add files to queue':{'cmd': self.addtoQueue}}       
        self.queue_menu=self.create_pulldown(self.menu,self.queue_menu)        
        self.menu.add_cascade(label='Queue',menu=self.queue_menu['var'])        
        self.help_menu={ '01Online Help':{'cmd': self.help} }       
        self.help_menu=self.create_pulldown(self.menu,self.help_menu)
        self.menu.add_cascade(label='Help',menu=self.help_menu['var'])        
        self.main.config(menu=self.menu)        
        return

    def presetsMenu(self):
        """Create preset entries"""
        self.preset_menu=Menu(self.menu,tearoff=0)
        self.menu.add_cascade(label='Presets',menu=self.preset_menu)
        presets = ['J-810_cd_bywavelength', 'J-810_cd_bytemp', 'nmr_titration_sparky']
        for name in presets:
            def func(p):
                def new():
                    self.p.loadPreset(preset=p)
                return new
            self.preset_menu.add_command(label=name,command=func(name))
        return       
       
    def openRaw(self, filename=None):
        """Open a raw file, if more than one file we add them to queue"""
        if filename==None:           
            filename = self.openFilename('.csv')                    
        if not os.path.exists(filename): return
        #now we open the first file only
        lines = self.p.openRaw(filename)
        self.updateinfoPane()
        self.showRawFile(lines)
        self.showPreview()       
        return
        
    def createConfig(self):               
        filename = self.saveFilename()
        self.p.createConfig(filename)
        return
            
    def loadConfig(self, filename=None):
        if filename==None:
            filename=self.openFilename('.conf')
        if not filename: return
        self.p.parseConfig(filename)
        self.updateinfoPane()
        return        
        
    def editConfig(self):
        self.editFile(self.p.conffile)
        return
                
    def editFile(self, filename=None):
        """Edit a file"""
        if filename==None:
            filename = self.openFilename('.conf')
        if not filename:
            return
        from PEATDB.textFrame import textFrame
        tf = textFrame(parent=self,title=filename)
        tf.load_from_file(filename)
        return
   
    def showRawFile(self, lines):
        """Show raw file contents"""
        #if self.p.rowend>len(lines):
        #    self.p.rowend=len(lines)
        c=self.rawcontents
        c.delete("1.0",END)
        c.component('columnheader').delete("1.0",END)
        c.component('rowheader').delete("1.0",END)
        count=0
        for row in range(0, len(lines)):
            if type(lines[row]) is types.StringType:
                line = string.strip(lines[row])
            else:
                line = lines[row]         
            c.insert(END,'%s\n' %line)
            c.component('rowheader').insert(END, str(count)+'\n')
            count=count+1
            
        return

    def execute(self):
        self.log.delete(1.0,END)
        self.p.run()
        return
    
    def showPreview(self,lines=None):
        """Show how the data looks with the import formatting applied"""
        self.previewer.update()        
        return        
        
    def addtoQueue(self, files=None):
        """Add files"""
        if files==None:
            files = self.openFilenames()
        self.p.addtoQueue(files)   
        self.updateinfoPane()
        self.queueFrame.update()
        return
    
    def openEkin(self, fname=None):
        """Open results in ekin"""
        
        if len(self.p.results)==0:
            print 'no results files'
            return
        fname = self.p.results[0]        
        EK = EkinApp(parent=self, project=fname)
        return  
    
    def openFilename(self, ext='.txt'):
        filename=tkFileDialog.askopenfilename(defaultextension=ext,initialdir=self.p.savedir,
                                              filetypes=[("text files","*.txt"),
                                                         ("csv files","*.csv"),
                                                         ("csvx files","*.csvx"),
                                                         ("excel files","*.xls"),
                                                         ("conf files","*.conf"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        return filename

    def openFilenames(self, ext='.txt'):
        filename=tkFileDialog.askopenfilenames(defaultextension=ext,initialdir=self.p.savedir,
                                              filetypes=[("text files","*.txt"),                                                                                                              
                                                         ("All files","*.*")],
                                              parent=self.main)
        return filename
    
    def saveFilename(self, ext='.conf'):
        filename=tkFileDialog.asksaveasfilename(defaultextension=ext,initialdir=self.p.savedir,
                                              filetypes=[("conf files","*.conf"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        return filename
    
    def write(self, txt):
        """Handle stdout"""
        self.log.appendtext(txt)       
        self.log.update_idletasks()
        return

    def help(self):
        import webbrowser
        link='http://enzyme.ucd.ie/main/index.php/EkinPipe'
        webbrowser.open(link,autoraise=1)        
        return
    
    def quit(self):
        self.main.destroy()
        if not self.parent:
            sys.exit()
        return

class PlotPreviewer(Frame):
    def __init__(self, parent=None, app=None):
        self.parent = parent 
        self.app = app
        self.p = self.app.p #reference to pipeline object
        if not self.parent:
            Frame.__init__(self)        
        
        fr = Frame(self)
        b=Button(fr,text='update',command=self.update)
        b.pack(side=TOP,fill=BOTH)
        b=Button(fr,text='prev',command=self.prev)
        b.pack(side=TOP,fill=BOTH)     
        b=Button(fr,text='next',command=self.next)
        b.pack(side=TOP,fill=BOTH)
        self.numplotscounter = Pmw.Counter(fr,
                labelpos='w',
                label_text='plots:',
                entryfield_value=1,
                entry_width=3,
                datatype = 'integer',
                entryfield_command=self.replot,
                entryfield_validate={'validator':'numeric', 'min':1,'max':8})    
        self.numplotscounter.pack()
        fr.pack(side=LEFT)
        self.plotframe = PlotPanel(parent=self, side=BOTTOM, height=200)
        self.dsindex = 0
        self.opts = {'markersize':15,'fontsize':9}
        return

    def replot(self):
        p = int(self.numplotscounter.getvalue())
        if p>1:
            dsets=self.E.datasets[self.dsindex:self.dsindex+p]                                      
            c=p/2
        else:    
            dsets = self.E.datasets[self.dsindex]
            c=1            
        self.plotframe.plotCurrent(dsets,options=self.opts,cols=c)       
    
    def loadData(self, data):
        """Load dict into datasets"""
        #table=self.previewTable = TableCanvas(frame, **self.tableformat)
        #table.createTableFrame()
        self.E = E = EkinProject(mode='General')
        for d in data.keys():
            xy = data[d]
            ek=EkinDataset(xy=xy)
            E.insertDataset(ek, d)
        self.plotframe.setProject(E)                   
        return
        
    def update(self, evt=None):        
        """Reload data dict from main app"""
        '''try:
            data = self.p.doImport()
        except:
            print 'could not do import with current config'
            data=None'''
        data = self.p.doImport()    
        if data == None: return
        self.dsindex = 0
        self.loadData(data)
        self.replot()
        return
        
    def prev(self):
        if self.dsindex <= 0:
            self.dsindex = 0
        else:    
            self.dsindex -= 1
        self.replot()    
        return   
        
    def next(self):
        if self.dsindex >= self.E.length-1:
            self.dsindex = self.E.length-1
        else:    
            self.dsindex += 1
        self.replot()        
        return
        
class queueManager(Frame):
    """Small class for managing file queue"""
    def __init__(self, parent=None, app=None):
        self.parent = parent 
        self.app = app
        self.p = self.app.p #reference to pipeline object
        if not self.parent:
            Frame.__init__(self)
        self.listbox = Pmw.ScrolledListBox(self,
                labelpos='n',
                label_text='File Queue',
                listbox_height = 8,                
                dblclickcommand=self.setFile,
                usehullsize = 1,
                hull_width = 200,
                hull_height = 150)
        self.listbox.component('listbox').configure(selectmode=EXTENDED)
        self.listbox.component('listbox').configure(selectbackground='lightblue')
        self.listbox.pack(fill=BOTH,expand=1)
        self.doButtons()
        return

    def doButtons(self):
        fr=Frame(self)
        methods = {'remove':self.removeSelected,'clear':self.clear}        
        for m in methods:
            b=Button(fr,text=m,command=methods[m])
            b.pack(side=LEFT,fill=BOTH)
        fr.pack(fill=BOTH,side=BOTTOM)    
        return
    
    def clear(self):
        self.listbox.delete(0,END)
        self.p.queue = []
        
    def removeSelected(self):
        items = self.listbox.getcurselection()      
        pos = 0
        for i in items:
            idx = items.index(i)
            self.listbox.delete( idx,idx )
            pos = pos + 1 
        return
        
    def update(self):
        print self.p.queue 
        flist=self.p.queue
        self.listbox.setlist(flist)
        
    def setFile(self):
        sel = self.listbox.getcurselection()[0]
        print sel
        self.app.openRaw(sel)
        
        
class Pipeline(object):    
    """This class does all the pipeline processing and configuration"""
    
    sections = ['base','fitting','excel']
    
    def __init__(self, conffile=None):
        
        if conffile==None:
            self.createConfig('pipe.conf')
        self.savedir = os.getcwd()
        self.filename = ''
        self.lines  =None
        self.queue = []
        self.results = []
        return
    
    def createConfig(self, filename, **kwargs):
        """Create a basic config file with default options"""
        
        c = ConfigParser.ConfigParser()
        s = 'base'
        c.add_section(s)
        c.set(s, 'format', 'databyrow')
        c.set(s, 'delimeter', ',')
        c.set(s, 'decimalsymbol', '.')
        c.set(s, 'checkunicode', 0)
        c.set(s, 'path', os.getcwd())
        c.set(s, 'columns', 10)
        c.set(s, 'rowstart', 0)
        c.set(s, 'colstart', 0)
        c.set(s, 'rowend', 1000)
        c.set(s, 'colnamesstart', 0)             
        c.set(s, 'groupbycol', 0)
        c.set(s, 'alternatecols', 0)
        c.set(s, 'colrepeat', 0)
        c.set(s, 'groupbyrow', 0)
        c.set(s, 'alternaterows', 0)
        c.set(s, 'rowrepeat', 0)
        c.set(s, 'rowdataformat', 'number')
        c.set(s, 'coldataformat', 'number')
        c.set(s, 'timeformat', "%Y-%m-%d %H:%M:%S")
        
        f = 'fitting'
        c.add_section(f)
        c.set(f, 'xerror', 0.1)
        c.set(f, 'yerror', 0.1)
        c.set(f, 'model1', 'linear')
        c.set(f, 'iterations1', 20)
        c.set(f, 'ignorecomments', 1)
        #excel settings
        e='excel'
        c.add_section(e)
        c.set(e, 'sheet', 0)
        c.set(e, 'numsheets', 1)
        #use kwargs to create specific settings
        for k in kwargs:
            c.set(s, k, kwargs[k])
        c.write(open(filename,'w'))
        self.parseConfig(filename)       
        return c

    def parseConfig(self, conffile=None):
        """Parse the config file"""        
        f = open(conffile,'r')
        cp = ConfigParser.ConfigParser()
        
        try:
            cp.read(conffile)
        except:
            print 'failed to read config file! check format'
            return
        self.conffile = conffile
        
        #get format
        fmt = cp.get('base', 'format')
        self.importer = self.getImporter(fmt, cp)
        print self.importer
        print 'parsed config file ok\n'
        return
         
    def getImporter(self, format, cp):
        """Get the required importer object"""
        if format == 'databyrow':
            importer = DatabyRowImporter(cp)
        elif format == 'databycolumn':
            importer = DatabyColImporter(cp)

        return importer
            
    def loadPreset(self, preset=None):
        """Load preset config for specific exp type data"""        
        if preset == 'J-810_cd_bywavelength':
            self.createConfig(preset+'.conf',format='databyrow',columns=10,rowstart=19,
                              delimeter='tab',xerror=0,yerror=0.2,
                              model1='Sigmoid')
        if preset == 'J-810_cd_bytemp':
            self.createConfig(preset+'.conf',format='databycolumn',columns=10,rowstart=19,rowend=421,
                              delimeter='tab',xerror=0,yerror=0.2,
                              datainrows=0,
                              model1='')            
        elif preset == 'nmr_titration_sparky':
            self.createConfig(preset+'.conf',columns=10,rowstart=0,
                              delimeter='tab',xerror=0.1,yerror=0.02,
                              model1='1 pKa 2 Chemical shifts')            
            
        print 'preset conf file written, you can also rename and edit this'
        return
        
    def openRaw(self, filename=None, callback=None):
        """Open raw file, display preview and get some info about them"""
        
        if os.path.splitext(filename)[1] == '.xls':
            lines = self.openExcel(filename)            
        else:
            fd=open(filename)
            lines = fd.readlines()
            fd.close()
        if lines == None:
            return None
        
        print 'opened file %s' %filename
        self.filename = filename
        self.queue = [filename]        
        self.filename = filename
        self.lines = lines
        return lines

    def openExcel(self, filename):
        """Open raw excel file"""
        
        try:
            import xlrd
        except:
            print 'xlrd required for excel import'
            return
        lines=[]
        sep=' '
        if filename != None:
            book = xlrd.open_workbook(filename)
            print 'The number of worksheets is', book.nsheets
            sh = book.sheet_by_index(self.sheet)           
            print 'current sheet: %s rows: %s cols: %s' %(sh.name, sh.nrows, sh.ncols)
            for r in range(sh.nrows):
                vals=sh.row_values(r)
                txt=''
                for v in vals:
                    txt=txt+str(v)+sep                
                lines.append(txt)
        return lines
        
    def doImport(self, lines=None):
        """Import file with current setting and return a dict"""
        if lines == None:
            if self.lines!=None:
                lines=self.lines
            else:
                print 'no file loaded yet'
                return None
        if self.conffile == None:
            self.loadConfig()
        else:    
            self.parseConfig(self.conffile)
        try:
            data = self.importer.doImport(lines)
        except Exception, e:
            print e
            data = None
        return data
        
    def run(self):
        """Do pipeline with the current config"""
        
        #clear log
        self.results = [] #list of files
        for filename in self.queue:
            lines = self.openRaw(filename)                      
            rawdict = self.doImport(lines)
            print 'processing raw data..'
            E = EkinProject(mode='General')
            for d in rawdict.keys():          
                xy = rawdict[d]
                ek=EkinDataset(xy=xy)
                E.insertDataset(ek, d)
            '''if self.model1 != '':    
                E.fitDatasets('ALL', models=[self.model1], noiter=self.iterations1, 
                               conv=1e-6, grad=1e-6, silent=True)
                for d in E.datasets:
                    ferrs = E.estimateExpUncertainty(d, runs=10)
                    E.addMeta(d, 'exp_errors', ferrs)'''
            prjname=self.filename+'.ekinprj'
            E.saveProject(prjname)
            self.results.append(prjname)
            print E, 'saved to %s' %prjname
        print 'done'
        return   
    
    def addtoQueue(self, files):
        """Add files"""     
     
        for f in files[:]:
            if f not in self.queue:
                self.queue.append(f)       
        return
            
class BaseImporter(object):
    """Importer class, sub-class this to define methods specific to each kind of
       import format"""
       
    def __init__(self, cp):
        """Arguments: 
            cp - a ConfigParser object that has been loaded in the parent app""" 
        for s in Pipeline.sections:
            for f in cp.items(s):
                #print f[0], f[1]
                try: val=int(f[1])
                except: val=f[1]
                self.__dict__[f[0]] = val
        if self.delimeter=='': self.delimeter=' '
        elif self.delimeter=='tab': self.delimeter='\t'                 
        return
           
    def getRow(self, row, lines):
        line = string.strip(lines[row]).split(self.delimeter)
        return row        
       
    def getColumn(self, c, lines):
        return

    def getRowHeader(self, lines):
        """Return labels in header column"""
        labels=[]
        for row in range(self.rowstart+1, self.rowend):
            #if row>=len(lines):
            #    break
            line = string.strip(lines[row]).split(self.delimeter)
            labels.append(line[self.colstart])
        return labels    
        
    def getColumnHeader(self, lines):
        """Return labels in row header"""
        labels = string.strip(lines[self.rowstart]).split(self.delimeter)
        #column names might be offset from column data
        if self.colnamesstart != 0:
            labels = labels[self.colnamesstart:]
        return labels    
         
    def checkValue(self, val):
        """Coerce a string to float if possible"""
        #add code to handle commas in thousand separators
        dec = self.decimalsymbol
        if dec == '.':
            try:
                return float(val)
            except:
                return None
        else:        
            try:
                return float(val.replace(".","").replace(dec,"."))
            except ValueError:
                return None
         
    def checkTime(self, val, timeformat):
        """Coerce to a datetime object"""
        try:
            datetime.strptime(val,timeformat)
        except:    
            return None
        
    def checkUnicode(self, s):
        """Check for unicode string"""
        try:
            s.decode('ascii')
        except UnicodeDecodeError:
            s = unicode(s)
        return s
        
    def doImport(self, lines):
        """Should be overrriden"""        
        return
 
class DatabyColImporter(BaseImporter):
    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):
        """Import where data are in cols"""
    
        data = {}        
        names = self.getColumnHeader(lines)

        #handle if data is grouped in rows  
        if self.groupbyrow == True:
            pass
        #handle case if rows are grouped in alternating rows

        print names
        for col in range(self.colstart+1, len(names)):
            name=names[col]
            x=[]; y=[]
            if self.alternaterows == True:
                step = self.rowrepeat
            else:
                step = 1
            
            for row in range(self.rowstart+1, self.rowend, step):  
                line = string.strip(lines[row]).split(self.delimeter)
                line = line[len(line)-len(names):]
                #print line, row
                if col >= len(line): continue
                a = self.checkValue(line[0])
                b = self.checkValue(line[col])
                #print a, b
                if a==None or b==None:
                    continue
                x.append(a); y.append(b)
            if len(x)<1: continue

            #print name, x
            data[name] = [x,y]        
        return data
        
class DatabyRowImporter(BaseImporter):
    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):
        """Import where data are in rows"""
        data = {}
        colnames = string.strip(lines[self.rowstart]).split(self.delimeter)
        for row in range(self.rowstart+1, self.rowend):
            if row>=len(lines):
                break
            line = string.strip(lines[row]).split(self.delimeter)
            name=line[0]
            if self.ignorecomments==True and name.startswith('#'):
                continue
            x=[]; y=[]
            for c in range(1,self.columns):
                col=float(colnames[c])
                x.append(col)
                y.append(float(line[c]))
            data[name] = [x,y]        
        return data
        
        
def main():  
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-c", "--conf", dest="conf",
                            help="Provide a conf file", metavar="FILE")
    parser.add_option("-f", "--file", dest="file",
                        help="Raw file", metavar="FILE")    
    
    opts, remainder = parser.parse_args()
    app = PipeApp(rawfile=opts.file)
    if opts.conf != None:
        app.loadConfig(opts.conf)    
    if opts.file != None:
        app.openRaw(opts.file)
    app.mainloop()    

if __name__ == '__main__':
    main()
