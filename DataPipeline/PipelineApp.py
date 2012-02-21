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

from Base import Pipeline
from PEATDB.Ekin.Base import *
from PEATDB.Tables import TableCanvas
from PEATDB.TableModels import TableModel
from PEATDB.Ekin.Ekin_main import EkinApp, PlotPanel
from PEATDB.Ekin.Pylab import Options
from PEATDB.Ekin.ModelDesign import ModelDesignApp
import os, sys, math, random, glob, numpy, string, types
import csv
from Tkinter import *
import Pmw
import tkFileDialog
from PEATDB.GUI_helper import *
from Helper import HelperDialog

class PipeApp(Frame, GUI_help):
    """Data pipe GUI for importing and fitting of raw data.
       This class uses ekin provided an automated pipeline for fitting
       raw text data files and propagating errors.
       Uses a config file to store the pipeline settings"""

    def __init__(self, parent=None, rawfile=None, conffile=None):

        self.parent=parent
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
        else:
            self.main=Toplevel()
            self.master=self.main
        self.main.title('DataPipeline Desktop')
        ws = self.main.winfo_screenwidth()
        hs = self.main.winfo_screenheight()
        w = 800; h=600
        x = (ws/2)-(w/2); y = (hs/2)-(h/2)
        self.main.geometry('%dx%d+%d+%d' % (w,h,x,y))
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
        self.infopane = Frame(self.main,height=20)
        self.infopane.pack(side=BOTTOM,fill=BOTH,pady=4)
        self.updateinfoPane()
        Label(self.infopane,text='Conf file:').pack(side=LEFT)
        Label(self.infopane,textvariable=self.conffilevar,fg='darkblue').pack(side=LEFT,padx=4)
        Label(self.infopane,text='Files in queue:').pack(side=LEFT,padx=4)
        Label(self.infopane,textvariable=self.queuefilesvar,fg='darkblue').pack(side=LEFT)
        Label(self.infopane,text='Current file:').pack(side=LEFT,padx=4)
        Label(self.infopane,textvariable=self.currfilevar,fg='darkblue').pack(side=LEFT)
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
        self.previewer = PlotPreviewer(self.m1,app=self)
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
        return

    def updateinfoPane(self):
        if hasattr(self.p, 'conffile'):
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
                        '02Run Config Helper':{'cmd': self.launchHelper},
                        '03Model Design':{'cmd': self.launchModelDesigner},
                        '04Launch Ekin':{'cmd':self.openEkin}}
        self.run_menu=self.create_pulldown(self.menu,self.run_menu)
        self.menu.add_cascade(label='Run',menu=self.run_menu['var'])
        self.queue_menu={'01Add files to queue':{'cmd': self.addtoQueue},
                         '02Add folder to queue':{'cmd': self.addFolder}}
        self.queue_menu=self.create_pulldown(self.menu,self.queue_menu)
        self.menu.add_cascade(label='Queue',menu=self.queue_menu['var'])
        self.utils_menu={'01Run Tests':{'cmd': self.runTests}}
        self.utils_menu=self.create_pulldown(self.menu,self.utils_menu)
        self.menu.add_cascade(label='Utilities',menu=self.utils_menu['var'])        
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
        if filename == None:
            filename = self.openFilename('.conf')
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
        """Run current files in queue"""
        if len(self.p.queue) == 0:
            return
        from Dialogs import ProgressDialog
        pb = ProgressDialog(self.main)
        #self.log.delete(1.0,END)
        self.p.run(callback=pb.bar.update)
        pb.close()
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

    def addFolder(self, path=None):
        if path==None:
            path = self.openDirectory()
        self.p.addFolder(path)
        self.updateinfoPane()
        self.queueFrame.update()
        return

    def runTests(self):
        """Run tests"""
        from Testing import Tester
        t=Tester()
        t.formatTests(t.basictests)
        print 'tests completed ok'
        return
        
    def openEkin(self, fname=None):
        """Open ekin"""
        EK = EkinApp(parent=self)
        return

    def launchModelDesigner(self):
        self.modelapp = ModelDesignApp(parent=self)
        return

    def launchHelper(self):
        wz = HelperDialog(parent=self)
        return

    def openFilename(self, ext='.txt'):
        filename=tkFileDialog.askopenfilename(defaultextension=ext,
                                              initialdir=self.p.savedir,
                                              filetypes=[("text files","*.txt"),
                                                         ("csv files","*.csv"),
                                                         ("csvx files","*.csvx"),
                                                         ("excel files","*.xls"),
                                                         ("conf files","*.conf"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        return filename

    def openFilenames(self, ext='.txt'):
        filename=tkFileDialog.askopenfilenames(defaultextension=ext,
                                              initialdir=self.p.savedir,
                                              filetypes=[("text files","*.txt"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        return filename

    def saveFilename(self, ext='.conf'):
        filename=tkFileDialog.asksaveasfilename(defaultextension=ext,
                                              initialdir=self.p.savedir,
                                              filetypes=[("conf files","*.conf"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        return filename

    def openDirectory(self):
        folder = tkFileDialog.askdirectory(parent=self.main,
                                            initialdir=os.getcwd(),
                                            title='Select folder')
        return folder

    def write(self, txt):
        """Handle stdout"""
        self.log.yview('moveto', '1')
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
    def __init__(self, master, app=None):

        Frame.__init__(self, master)
        self.app = app
        if hasattr(self.app,'p'):
            self.p = self.app.p     #reference to pipeline object
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
        E = self.E = Pipeline.getEkinProject(data)
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
        flist=self.p.queue
        self.listbox.setlist(flist)

    def setFile(self):
        sel = self.listbox.getcurselection()[0]
        print sel
        self.app.openRaw(sel)

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-c", "--conf", dest="conf",
                            help="Provide a conf file", metavar="FILE")
    parser.add_option("-f", "--file", dest="file",
                        help="Raw file", metavar="FILE")
    parser.add_option("-d", "--dir", dest="directory",
                        help="Folder of raw files")

    opts, remainder = parser.parse_args()
    app = PipeApp(rawfile=opts.file)
    if opts.conf != None:
        app.loadConfig(opts.conf)
    if opts.file != None:
        app.openRaw(opts.file)
    if opts.directory != None:
        app.addFolder(opts.directory)
    app.mainloop()

if __name__ == '__main__':
    main()
