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
from PEATDB.TableModels import TableModel
import os, sys, math, random, glob, numpy, string
import ConfigParser, csv
from Tkinter import *
import Pmw
import tkFileDialog
from PEATDB.GUI_helper import *

class EkinPipe(Frame, GUI_help):
    """Data pipe GUI for importing and fitting of raw data.
       This class uses ekin provided an automated pipeline for fitting
       raw text data files abnd propagating errors.
       Uses a config file to store the pipeline settings"""
    
    def __init__(self, parent=None, rawfile=None, conffile=None):
        self.parent=parent
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
        else:
            self.main=Toplevel()
            self.master=self.main
        self.main.title('EkinPipe - Fitting Pipeline')
        self.main.geometry('800x600+200+100')
        self.main.protocol('WM_DELETE_WINDOW',self.quit)
        self.savedir = os.getcwd()
        self.filename=''
        self.queue = []
        self.results = []
        self.setupVars()
        if conffile==None:
            self.createConfig('pipe.conf')
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
        self.m.pack(side=BOTTOM,fill=BOTH,expand=1) 
        self.preview = Pmw.ScrolledText(self.m,
                labelpos = 'n',
                label_text='Raw File Preview',
                rowheader=1,
                columnheader=1,
                Header_foreground = 'blue',
                rowheader_width = 3,
                usehullsize = 1,
                hull_width = 500,
                hull_height = 500,
                text_wrap='none')
        self.m.add(self.preview)
        self.log = Pmw.ScrolledText(self.m,
                labelpos = 'n',
                label_text='Logs',
                usehullsize = 1,
                hull_width = 400,
                hull_height = 500,
                text_wrap='word')       
        self.m.add(self.log)
        self.infopane = Frame(self.main,height=20)
        self.infopane.pack(side=BOTTOM,fill=X,pady=4)        
        self.updateinfoPane()
        Label(self.infopane,text='Conf file:').pack(side=LEFT)
        Label(self.infopane,textvariable=self.conffilevar,fg='darkblue').pack(side=LEFT,padx=4)
        Label(self.infopane,text='Files in queue:').pack(side=LEFT,padx=4)
        Label(self.infopane,textvariable=self.queuefilesvar,fg='darkblue').pack(side=LEFT)
        Label(self.infopane,text='Current file:').pack(side=LEFT,padx=4)
        Label(self.infopane,textvariable=self.currfilevar,fg='darkblue').pack(side=LEFT)        
        return

    def updateinfoPane(self):
        self.conffilevar.set(self.conffile)
        self.queuefilesvar.set(len(self.queue))
        self.currfilevar.set(self.filename)
        return
    
    def createMenuBar(self):
        """Create the menu bar for the application. """
        self.menu=Menu(self.main)
        self.file_menu={ '01Open Raw File(s)':{'cmd':self.openRaw},
                         '02Load Config File':{'cmd':self.loadConfig},
                         '03Edit Current Config':{'cmd':self.editConfig},
                         '04Create Config':{'cmd':self.createConfig},                        
                         '05Quit':{'cmd':self.quit}}                         
        self.file_menu=self.create_pulldown(self.menu,self.file_menu)
        self.menu.add_cascade(label='File',menu=self.file_menu['var'])
        self.presetsMenu()                
        self.run_menu={ '01Execute':{'cmd': self.run},
                        '02Add files to queue':{'cmd': self.addtoQueue},
                        '03Open results in Ekin':{'cmd':self.openEkin}}       
        self.run_menu=self.create_pulldown(self.menu,self.run_menu)
        self.menu.add_cascade(label='Run',menu=self.run_menu['var'])
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
                    self.loadPreset(preset=p)
                return new
            self.preset_menu.add_command(label=name,command=func(name))
        return
    
    def loadConfig(self):
        f=self.openFilename()
        if not f: return
        self.parseConfig(f)
        return
        
    def editConfig(self):
        self.editFile(self.conffile)
        return
        
    def parseConfig(self, conffile=None):
        """Parse a config file - settings for import/fitting"""        
        f = open(conffile,'r')
        c = ConfigParser.ConfigParser()
        try:
            c.read(conffile)
        except:
            print 'failed to read config file!'
            return
        self.conffile = conffile
        for f in c.items('settings'):
            print f[0], f[1]
            try: val=int(f[1])
            except: val=f[1]
            self.__dict__[f[0]] = val                
        print 'parsed config file ok\n'
        self.updateinfoPane()
        return
    
    def createConfig(self, filename=None, **kwargs):
        """Create a basic config file with default options"""        
        c = ConfigParser.ConfigParser()
        s = 'settings'
        c.add_section(s)
        c.set(s, 'delimeter', ',')
        c.set(s, 'columns', 10)
        c.set(s, 'rowstart', 0)
        c.set(s, 'colstart', 0)
        c.set(s, 'rowend', 1000)
        c.set(s, 'rowrep', 0)
        c.set(s, 'datainrows', 1)
        c.set(s, 'path', os.getcwd())        
        c.set(s, 'xerror', 0.1)
        c.set(s, 'yerror', 0.1)
        c.set(s, 'model1', 'linear')
        c.set(s, 'iterations1', 20)
        c.set(s, 'ignorecomments', 1)
        #excel settings
        c.set(s, 'sheet', 0)
        c.set(s, 'numsheets', 1)
        #use kwargs to create specific settings
        for k in kwargs:
            c.set(s, k, kwargs[k])
        #dump to file
        if filename == None:
            filename = self.saveFilename()
        c.write(open(filename,'w'))
        self.parseConfig(filename)       
        return c

    def loadPreset(self, preset=None):
        """Load preset config for specific exp type data"""        
        if preset == 'J-810_cd_bywavelength':
            self.createConfig(preset+'.conf',columns=10,rowstart=19,
                              delimeter='tab',xerror=0,yerror=0.2,
                              model1='Sigmoid')
        if preset == 'J-810_cd_bytemp':
            self.createConfig(preset+'.conf',columns=10,rowstart=19,rowend=421,
                              delimeter='tab',xerror=0,yerror=0.2,
                              datainrows=0,
                              model1='')            
        elif preset == 'nmr_titration_sparky':
            self.createConfig(preset+'.conf',columns=10,rowstart=0,
                              delimeter='tab',xerror=0.1,yerror=0.02,
                              model1='1 pKa 2 Chemical shifts')            
            
        print 'preset conf file written, you can also rename and edit this'
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

    def openRaw(self, filename=None):
        """Open raw file, display preview and get some info about them"""
        if filename==None:
            filename = self.openFilename('.csv')
        if not filename or not os.path.exists(filename):
            return   
        if os.path.splitext(filename)[1] == '.xls':
            lines = self.openExcel(filename)            
        else:
            fd=open(filename)
            lines = fd.readlines()
            fd.close()
        if lines == None:
            return None
        self.showPreview(lines)
        print 'opened file %s' %filename
        self.filename = filename
        self.queue = [filename]
        self.updateinfoPane()
        return lines

    def openExcel(self, filename):
        """Open raw excel file"""
        try:
            import xlrd
        except:
            print 'xlrd required for excel import'
            return
        lines=[]
        sep=','
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
    
    def showPreview(self, lines):
        """Update preview of raw file"""
        if self.rowend>len(lines):
            self.rowend=len(lines)
        p=self.preview
        p.delete("1.0",END)
        p.component('columnheader').delete("1.0",END)
        p.component('rowheader').delete("1.0",END)
        count=0
        for row in range(0, len(lines)):
            if type(lines[row]) is types.StringType:
                line = string.strip(lines[row])
            else:
                line = lines[row]         
            p.insert(END,'%s\n' %line)
            p.component('rowheader').insert(END, str(count)+'\n')
            count=count+1
        return

    def addtoQueue(self):
        """Add files"""        
        files = self.openFilenames()
        for f in files[:]:
            if f not in self.queue:
                self.queue.append(f)        
        print self.queue
        self.updateinfoPane()
        return

    def doImport(self, filename=None):
        """Import file with current setting and send to dict"""
        if self.conffile == None:
            self.loadConfig()
        else:    
            self.parseConfig(self.conffile)        
        data = {}
        lines = self.openRaw(filename)
        if self.delimeter=='': self.delimeter=' '
        elif self.delimeter=='tab': self.delimeter='\t'
        if self.datainrows == 1:
            #datasets per row            
            colnames = string.strip(lines[self.rowstart]).split(self.delimeter)   
            print colnames
            for row in range(self.rowstart+1, self.rowend):
                if row>=len(lines):
                    break
                line = string.strip(lines[row]).split(self.delimeter)
                name=line[0]
                if self.ignorecomments==1 and name.startswith('#'):
                    continue               
                x=[]; y=[]
                for c in range(1,self.columns):
                    col=float(colnames[c])
                    x.append(col)
                    y.append(float(line[c]))
                data[name] = [x,y]
        else:
            #datasets per column
            names = string.strip(lines[self.rowstart]).split(self.delimeter)
            colnames = ['x','y']
            for col in range(self.colstart+1, len(names)):
                name=names[col]
                x=[]; y=[]
                for row in range(self.rowstart+1, self.rowend):
                    if row>=len(lines):
                        break                    
                    line = string.strip(lines[row]).split(self.delimeter)
                    try:
                        x.append(float(line[0]))
                        y.append(float(line[col]))
                    except:
                        continue
                #print name, len(x),len(y)    
                data[name] = [x,y]
        print 'imported data ok, found %s datasets' %len(data)
        return data
    
    def run(self):
        """Do pipeline with current config"""
        #clear log
        self.results = [] #list of files
        self.log.delete(1.0,END)
        for filename in self.queue:
            raw = self.doImport(filename)
            print 'processing raw data..'
            E = EkinProject(mode='General')
            for d in raw.keys():          
                xy = raw[d]
                ek=EkinDataset(xy=xy)
                E.insertDataset(ek, d)
            if self.model1 != '':    
                E.fitDatasets('ALL', models=[self.model1], noiter=self.iterations1, conv=1e-6, grad=1e-6, silent=True)
                '''for d in E.datasets:
                    ferrs = E.estimateExpUncertainty(d, runs=10)
                    E.addMeta(d, 'exp_errors', ferrs)'''
            prjname=self.filename+'.ekinprj'
            E.saveProject(prjname)
            self.results.append(prjname)
            print E, 'saved to %s' %prjname
        print 'done'
        return

    def openEkin(self, fname=None):
        """Open results in ekin"""
        
        if len(self.results)==0:
            print 'no results files'
            return
        fname = self.results[0]
        from Ekin_main import EkinApp
        EK = EkinApp(parent=self, project=fname)
        return
    
    def openFilename(self, ext='.txt'):
        filename=tkFileDialog.askopenfilename(defaultextension=ext,initialdir=self.savedir,
                                              filetypes=[("text files","*.txt"),
                                                         ("csv files","*.csv"),
                                                         ("csvx files","*.csvx"),
                                                         ("excel files","*.xls"),
                                                         ("conf files","*.conf"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        return filename

    def openFilenames(self, ext='.txt'):
        filename=tkFileDialog.askopenfilenames(defaultextension=ext,initialdir=self.savedir,
                                              filetypes=[("text files","*.txt"),                                                                                                              
                                                         ("All files","*.*")],
                                              parent=self.main)
        return filename
    
    def saveFilename(self, ext='.conf'):
        filename=tkFileDialog.asksaveasfilename(defaultextension=ext,initialdir=self.savedir,
                                              filetypes=[("conf files","*.conf"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        return filename
    
    def write(self, txt):
        """Handle stdout"""
        self.log.insert(END, txt)
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
        
def test():
    """Do basic tests"""

    return

def main():  
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-c", "--conf", dest="conf",
                            help="Provide a conf file", metavar="FILE")
    parser.add_option("-f", "--file", dest="file",
                        help="Raw file", metavar="FILE")    
    
    opts, remainder = parser.parse_args()
    app = EkinPipe(rawfile=opts.file)
    if opts.conf != None:
        app.parseConfig(opts.conf)    
    if opts.file != None:
        app.openRaw(opts.file)
    app.mainloop()    

if __name__ == '__main__':
    main()
