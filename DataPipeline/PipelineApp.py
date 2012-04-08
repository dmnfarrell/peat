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
from PEATDB.Ekin.Ekin_main import EkinApp
from PEATDB.Ekin.Plotting import PlotPanel, Options
from PEATDB.Ekin.ModelDesign import ModelDesignApp
import os, sys, math, random, glob, numpy, string, types
import pickle
from Tkinter import *
import Pmw
import tkFileDialog
from PEATDB.GUI_helper import *
from Helper import HelperDialog
from Editor import TextEditor
from Rename import BatchRenameApp
import Images

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
            self.conffilevar.set(self.p.configurationfile)
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
        self.project_menu={'01Load Project':{'cmd': self.loadProject},
                            '02Save Project':{'cmd': self.saveProject}}
        self.project_menu=self.create_pulldown(self.menu,self.project_menu)
        self.menu.add_cascade(label='Project',menu=self.project_menu['var'])
        self.run_menu={'01Execute':{'cmd': self.execute}}
        self.run_menu=self.create_pulldown(self.menu,self.run_menu)
        self.menu.add_cascade(label='Run',menu=self.run_menu['var'])
        self.queue_menu={'01Add files to queue':{'cmd': self.addtoQueue},
                          '02Add folder to queue':{'cmd': self.addFolder},
                          '03Clear queue':{'cmd': self.clearQueue}}
        self.queue_menu=self.create_pulldown(self.menu,self.queue_menu)
        self.menu.add_cascade(label='Queue',menu=self.queue_menu['var'])
        self.utils_menu={'01Show Config Helper':{'cmd': self.launchHelper},
                         '02Model Design':{'cmd': self.launchModelDesigner},
                         '03Launch Ekin':{'cmd':self.openEkin},
                         '04Text Editor':{'cmd': self.startTextEditor},
                         '05Batch File Rename':{'cmd': self.batchFileRename},
                         '06Clear Log':{'cmd': self.clearLog},
                         '07Run Tests':{'cmd': self.runTests}}
        self.utils_menu=self.create_pulldown(self.menu,self.utils_menu)
        self.menu.add_cascade(label='Utilities',menu=self.utils_menu['var'])
        self.help_menu={ '01Online Help':{'cmd': self.help},
                          '02About':{'cmd': self.about},}
        self.help_menu=self.create_pulldown(self.menu,self.help_menu)
        self.menu.add_cascade(label='Help',menu=self.help_menu['var'])
        self.main.config(menu=self.menu)
        return

    def openRaw(self, filename=None):
        """Open a raw file, if more than one file we add them to queue"""
        if filename==None:
            filename = self.openFilename()
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
            filename = self.openFilename('conf')
        if not filename: return
        self.p.parseConfig(filename)
        self.updateinfoPane()
        return

    def editConfig(self):
        self.editFile(self.p.configurationfile)
        return

    def editFile(self, filename=None):
        """Edit a file"""
        if filename==None:
            filename = self.openFilename('conf')
        if not filename:
            return

        tf = TextEditor(parent=self,title=filename)
        tf.openFile(filename)
        return

    def showRawFile(self, lines):
        """Show raw file contents"""

        c = self.rawcontents
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

    def stopCurrent(self):
        self.p.stop = True
        print 'cancel pressed.. please wait'
        return

    def execute(self):
        """Run current files in queue"""
        if len(self.p.queue) == 0:
            return
        from Dialogs import ProgressDialog
        signal=True

        self.pb = ProgressDialog(self.main, cancel=self.stopCurrent)
        self.pb.after(100, self.pb.updateValue())
        self.p.run(callback=self.pb.updateValue)
        self.pb.close()
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

    def updateFromQueue(self):
        """Check current file open after a queue deletion"""
        if self.p.filename not in self.p.queue.values():
            self.clearCurrentFile()
        return

    def clearCurrentFile(self):
        """Clear current file"""
        self.p.closeFile()
        self.updateinfoPane()
        self.rawcontents.clear()
        self.previewer.clear()
        return

    def addFolder(self, path=None):
        if path==None:
            path = self.openDirectory()
        self.p.addFolder(path)
        self.updateinfoPane()
        self.queueFrame.update()
        return

    def clearQueue(self):
        self.queueFrame.clear()
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
        if self.p.modelsfile != '':
            self.modelapp.loadModelsFile(self.p.modelsfile)
        return

    def launchHelper(self):
        wz = HelperDialog(parent=self)
        return

    def openFilename(self, ext=['txt','csv','xls']):
        if not type(ext) is types.ListType:
            ext=[ext]
        filetypes=[]
        for e in ext:
            filetypes.append(('%s files' %e,'*.%s' %e))
        filetypes.append(("All files","*.*"))
        filename=tkFileDialog.askopenfilename(defaultextension=ext,
                                              initialdir=self.p.savedir,
                                              filetypes=filetypes,
                                              parent=self.main)
        return filename

    def openFilenames(self, ext='txt'):
        filetypes = [('%s files' %ext,'*.%s' %ext)]
        filetypes.append(("All files","*.*"))
        filename=tkFileDialog.askopenfilenames(defaultextension=ext,
                                              initialdir=self.p.savedir,
                                              filetypes=filetypes,
                                              parent=self.main)
        return filename

    def saveFilename(self, ext=''):
        if ext!='':
            filetypes = [('%s files' %ext,'*.%s' %ext)]
        else:
            filetypes = []
        filetypes.append(("All files","*.*"))
        filename=tkFileDialog.asksaveasfilename(defaultextension='.'+ext,
                                              initialdir=self.p.savedir,
                                              filetypes=filetypes,
                                              parent=self.main)
        return filename

    def openDirectory(self):
        folder = tkFileDialog.askdirectory(parent=self.main,
                                            initialdir=os.getcwd(),
                                            title='Select folder')
        return folder

    def loadProject(self, filename=None):
        if filename == None:
            filename = self.openFilename('proj')
        f = open(filename,'r')
        try:
            self.p = pickle.load(f)
        except Exception,e:
            print 'failed to load project'
            print 'Error returned:', e
            return
        print 'loaded project', filename
        name = os.path.splitext(filename)[1]
        self.p.writeConfig(filename='%s.conf' %name)
        self.updateinfoPane()
        self.queueFrame.update()
        self.previewer.update()
        if self.p.lines != None:
            self.showRawFile(self.p.lines)
        return

    def saveProject(self, filename=None):
        """Save project file"""
        if filename == None:
            filename = self.saveFilename('proj')
        f = open(filename,'w')
        pickle.dump(self.p,f)
        f.close()
        return

    def startTextEditor(self):
        t = TextEditor(parent=self)
        return

    def batchFileRename(self):
        B = BatchRenameApp(parent=self)
        return

    def write(self, txt):
        """Handle stdout"""
        self.log.yview('moveto', '1')
        self.log.appendtext(txt)
        self.log.update_idletasks()
        #update the progress bar aswell
        if hasattr(self, 'pb'):
            self.pb.updateValue()
        return

    def clearLog(self):
        self.log.delete(1.0,END)
        return

    def help(self):
        import webbrowser
        link='http://enzyme.ucd.ie/main/index.php/DataPipeline'
        webbrowser.open(link,autoraise=1)
        return

    def about(self):
        win=Toplevel()
        win.geometry('+500+350')
        win.title('About DataPipeline')
        logo = Images.logo()
        label = Label(win,image=logo)
        label.image = logo
        label.pack(fill=BOTH,padx=4,pady=4)
        text="""DataPipeline App, Version 1.0
             is a python desktop and web application
             that uses a configuration file to automate the import of
             raw data in a variety of formats.\n
             Released under GPL v3\n
             (C) Copyright 2011- Damien Farrell """
        text=text.replace('\t','')
        text= ' '.join(text.split())
        Label(win,text=text,wraplength=400).pack(fill=Y,side=TOP,pady=4)
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
        self.E = None
        if hasattr(self.app,'p'):
            self.p = self.app.p     #reference to pipeline object
        fr = Frame(self)
        b=Button(fr,text='update',command=self.update)
        b.pack(side=TOP,fill=BOTH)
        self.previmg = Images.prev()
        b=Button(fr,text='prev',image=self.previmg,compound='left',command=self.prev)
        b.pack(side=TOP,fill=BOTH)
        self.nextimg = Images.next()
        b=Button(fr,text='next',image=self.nextimg,compound='left',command=self.next)
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
        self.overlayvar = BooleanVar(); self.overlayvar.set(False)
        Checkbutton(fr, text='overlay plots', variable=self.overlayvar, command=self.replot).pack(side=TOP,fill=BOTH)
        fr.pack(side=LEFT)
        self.plotframe = PlotPanel(parent=self, side=BOTTOM, height=200, tools=True)
        self.dsindex = 0
        self.plotframe.Opts.opts['fontsize']=10

        b=Button(fr,text='open in Ekin',command=self.loadEkin)
        b.pack(side=TOP,fill=BOTH)
        return

    def replot(self):
        """Replot"""
        p = int(self.numplotscounter.getvalue())
        if p>1:
            dsets = self.E.datasets[self.dsindex:self.dsindex+p]
            c=p/2            
        else:
            dsets = self.E.datasets[self.dsindex]
            c=1

        if self.overlayvar.get() == True:
            plotopt = 3
            self.plotframe.Opts.opts['title']=' '
            self.plotframe.Opts.opts['legend'] = 1
        else:
            plotopt = 2
            self.plotframe.Opts.opts['title']=None
        self.plotframe.plotCurrent(dsets,
                                    cols=c, plotoption=plotopt)
        return

    def loadData(self, data):
        """Load dict into datasets"""
        E = self.E = Pipeline.getEkinProject(data)
        self.plotframe.setProject(E)
        return

    def update(self, evt=None):
        """Reload data dict from main app"""
        self.p = self.app.p
        data = self.p.doImport()
        if data == None: return
        self.dsindex = 0
        self.loadData(data)
        self.replot()
        return

    def clear(self):
        self.E = None
        self.plotframe.clear()
        return

    def prev(self):
        if self.E == None: return
        if self.dsindex <= 0:
            self.dsindex = 0
        else:
            self.dsindex -= 1
        self.replot()
        return

    def next(self):
        if self.E == None: return
        if self.dsindex >= self.E.length-1:
            self.dsindex = self.E.length-1
        else:
            self.dsindex += 1
        self.replot()
        return

    def loadEkin(self):
        if self.E != None:
            EK = EkinApp(parent=self, project=self.E)
        return

class queueManager(Frame):
    """Small class for managing file queue"""
    def __init__(self, parent=None, app=None):
        self.parent = parent
        self.app = app
        #self.p = self.app.p #reference to pipeline object
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
        p = self.app.p
        self.listbox.delete(0,END)
        p.queue = []
        return

    def removeSelected(self):
        """Remove selected items from queue"""
        p = self.app.p
        if len(self.p.queue) == 0:
            return
        p.queue.values()
        selected = self.listbox.getcurselection()
        for key, value in p.queue.items():
            if value in selected:
                del self.p.queue[key]
        self.update()
        self.app.updateFromQueue()
        return

    def update(self):
        """Update queue listbox"""
        p = self.app.p
        flist = p.queue.values()
        self.listbox.setlist(flist)
        return

    def setFile(self):
        sel = self.listbox.getcurselection()[0]
        self.app.openRaw(sel)
        return


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-c", "--conf", dest="conf",
                            help="Provide a conf file", metavar="FILE")
    parser.add_option("-f", "--file", dest="file",
                        help="Raw file", metavar="FILE")
    parser.add_option("-d", "--dir", dest="directory",
                        help="Folder of raw files")
    parser.add_option("-p", "--project", dest="project",
                        help="Project file", metavar="FILE")

    opts, remainder = parser.parse_args()
    app = PipeApp(rawfile=opts.file)
    if opts.conf != None:
        app.loadConfig(opts.conf)
    if opts.file != None:
        app.openRaw(opts.file)
    if opts.directory != None:
        app.addFolder(opts.directory)
    if opts.project != None:
        app.loadProject(opts.project)
    app.mainloop()

if __name__ == '__main__':
    main()
