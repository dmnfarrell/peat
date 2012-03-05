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

try:
    from Plugins import Plugin
except:
    from PEATDB.Plugins import Plugin
import os, types, copy, pickle
from Tkinter import *
import Pmw
import PEATSA.WebApp.Data
import PEATSA.WebApp.UtilityFunctions
import PEATSA.Core as Core
from PEATDB.Dialogs import MultipleValDialog
from PEATDB.Actions import DBActions
from PEATDB.TableModels import TableModel
from PEATDB.Tables import TableCanvas
import tkMessageBox, tkSimpleDialog, tkFileDialog

class PEATSAPlugin(Plugin):
    """Template GUI plugin for PEAT App"""
    capabilities = ['gui']
    requires = ['PEATSA']
    menuentry = 'PEATSA Plugin'
    gui_methods = {'fetchJob':'Fetch Job from Server',                  
                   'editConfigFile' : 'Configure Server',
                   'help':'Help',
                   'quit':'Close Window'}
    buttonorder = ['createJobDialog','fetchJob','editConfigFile','help','quit']
    about = 'This plugin allows you to call PEATSA'

    calctypes = ['stability','binding','pka']
    
    def main(self, parent=None, DB=None):
      
        if parent == None:
            if DB != None:
                self.DB = DB
                self.setupConnection()
            else:
                return
        else:
            self.parent = parent
            self.DB = parent.DB
            if self.DB == None:
                self.displayNoDBWarning()
                return
            self._doFrame()
            self.setupConnection()
            print 'Updating jobs table..'
            self.updateJobs()
        return self

    def setupConnection(self):
        """Set up connection"""
        homepath = os.path.expanduser("~") 
        self.confpath = os.path.join(homepath, 'peatsa.conf')
        if os.path.exists(self.confpath):    
            configuration = Core.Environment.Configuration(filename=self.confpath)
        else:
            configuration = Core.Environment.Configuration(searchDefaultLocations=False)
            configuration.add_section('DATABASE')
            configuration.set('DATABASE', 'database', 'DBSAInterface')
            configuration.set('DATABASE', 'host', 'enzyme.ucd.ie')
            configuration.set('DATABASE', 'user', 'peatdb')
            configuration.set('DATABASE', 'password', '123')
            configuration.writeToFile(self.confpath)
            if self.parent != None:               
                tkMessageBox.showwarning("Connection Error",
                         'No PEATSA server configured, press configure server'
                         ' to set a server, username and password.')
        self.connect(configuration)
        return
    
    def _doFrame(self):
        self.mainwin = self.parent.createChildFrame(width=460,title='PEATSA Plugin')
        #self.mainwin = self.parent.create
        methods = self._getmethods()
        methods = [m for m in methods if m[0] in self.gui_methods.keys()]        
        l=Label(self.mainwin, text='PEATSA Interface')
        l.pack(side=TOP,fill=BOTH)
        self.tf=LabelFrame(self.mainwin,text='Project Calculations')
        self.tf.pack(side=TOP,fill=BOTH,expand=1)
        self.manageJobsButtons(self.mainwin)
        self._createButtons(methods)
        self.log = self.createLogWin(self.mainwin)
        self.log.pack(side=TOP,fill=BOTH,expand=1)
        self.stdout2Log()
        self.mainwin.bind("<Destroy>", self.quit)
        #self.parent.sidepane.bind("<Destroy>", self.test1)
        return

    def _createButtons(self, methods):
        """Dynamically create buttons for supplied methods, which is a tuple
            of (method name, label)"""
        mbutton=Menubutton(self.mainwin, text='Options', width=12,
                                         borderwidth=2, relief=RIDGE,
                                         activeforeground='red')
        menu=Menu(mbutton,tearoff=0)
        mbutton['menu']=menu
        mbutton.pack(side=BOTTOM,fill=BOTH)
        for m in methods:
            menu.add_radiobutton(label=self.gui_methods[m[0]], 
                                        indicatoron=0,                                       
                                        command=m[1])
        b=Button(self.mainwin,text='Create Calculation',command=self.createJobDialog)
        b.pack(side=BOTTOM,fill=BOTH)            
        return

    def updateJobsTable(self):
        """Show table for current jobs list"""
        self.checkJobsDict()
        jobdict = self.DB.meta.peatsa_jobs      
        M = TableModel()
        #open job log from file
        f=open('jobstates.log','r')
        jl = pickle.load(f) 
        for j in jobdict:            
            jobid = jobdict[j]           
            try:
                M.addRecord(j,state=jl[jobid]['State'],date=jl[jobid]['Date'])
            except:
                M.addRecord(j,state='Not in DB')
        self.jobstable = TableCanvas(self.tf, model=M, height=100, editable=False)
        self.jobstable.createTableFrame()       
        self.log.yview('moveto', 1)
        f.close()
        return

    def manageJobsButtons(self, parent):
        fr1 = Frame(parent)
        Button(fr1,text='View Results',command=self.showAllResults,bg='#ccFFFF').pack(side=TOP,fill=BOTH,expand=1)
        fr1.pack(fill=BOTH)
        Button(fr1,text='Merge Results',command=self.mergeCurrent).pack(side=TOP,fill=BOTH,expand=1)
        fr1.pack(fill=BOTH)        
        fr = Frame(parent)
        c='#ADD8E6'
        Button(fr,text='Show Details',command=self.viewDetails,bg=c).pack(side=LEFT,fill=BOTH,expand=1)
        Button(fr,text='Manage Results',command=self.manageResults,bg=c).pack(side=LEFT,fill=BOTH,expand=1)
        Button(fr,text='Remove',command=self.removeJob,bg=c).pack(side=LEFT,fill=BOTH,expand=1)
        Button(fr,text='Resubmit',command=self.resubmitJob,bg=c).pack(side=LEFT,fill=BOTH,expand=1)
        fr.pack(fill=BOTH)
        return 
        
    def createLogWin(self, parent):
        log = Pmw.ScrolledText(parent,
                borderframe=1,
                labelpos = 'n',
                label_text='Log',
                usehullsize = 1,
                hull_width = 800,
                hull_height = 200,
                text_wrap='word')        
        return log

    def stdout2Log(self):
        """Redirect stdout to app control"""
        sys.stdout = self
        sys.stderr = self
        return
       
    def log2Stdout(self):
        """return to stdout"""
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        return
 
    def write(self, txt):
        """Handle stdout if required"""        
        self.log.appendtext(txt)
        self.log.update_idletasks()
        return

    def flush(self):
        return
 
    def connect(self, configuration):
        """Create connection"""
        self.connection = PEATSA.WebApp.UtilityFunctions.ConnectionFromConfiguration(configuration)
        self.jobManager = PEATSA.WebApp.Data.JobManager(self.connection)
        self.jobManager.setJobStateLogging('jobstates.log',interval=60)
        print '\nConnection to server made sucessfully.\n'     
        return
        
    def createMutationList(self, filename=None):        
        self.mutationList = Core.Data.MutationListFile(create=False)        
        return

    def fetchJob(self):
        """Get job from it's db ID and add to list"""
        
        mpDlg = MultipleValDialog(title='Get Job',
                                    initialvalues=('','my job1'),
                                    labels=('ID','Your label',),
                                    types=('string','string'),
                                    parent=self.mainwin)
        if mpDlg.result == True:
            jobid = mpDlg.results[0]
            name = mpDlg.results[1]
        else:
            return
        job = PEATSA.WebApp.Data.Job(jobid, self.connection)       
        if job != None: 
            print 'adding job id %s to list' %job.identification
            self.storeJob(name, job)
            self.updateJobs()
        return

    def writetempPDB(self,name=None,pdbfile='refprot.pdb'):
        if name==None:
            name = self.DB.meta.refprotein
        pdblines = self.DB[name].Structure
        #pdbfile = 'refprot.pdb'
        fd=open(pdbfile,'w')
        for line in pdblines:
            fd.write(line)
        fd.close()
        return pdbfile
            
    def getrefPDBName(self):
        name = self.DB.meta.refprotein
        if self.DB[name].has_key('pdbname'):
            name = self.DB[name]['pdbname']
            return name.split('.')[0]            
        else:
            return ''
        
    def createJobDialog(self):
        """Get details from user using dialog
           required: structure, mutations, type of calc and a tag (optional)"""

        def validatename(text):
            if not hasattr(self.DB.meta,'peatsa_jobs'):
                return 1
            if text in self.DB.meta.peatsa_jobs:
              return -1
            else:
              return 1
        def close():
            jobdlg.destroy()
        def loadmuts():            
            filename=tkFileDialog.askopenfilename(initialdir=os.getcwd(),
                                       filetypes=[("All files","*")])
            if filename:
                mutlist.importfile(filename)
            return
        def loadmutsfromDB():
            for p in self.DB.getRecs():               
                mut = self.DB.get(p).Mutations
                if mut == None or mut=='':
                    continue
                if type(mut) is types.StringType:
                    mutlist.appendtext(mut+'\n')
                else:
                    mutstring = mut.getMutationString()
                    if mutstring != None:
                        mutlist.appendtext(mutstring+'\n')
            return
        def getstruct():
            filename=tkFileDialog.askopenfilename(defaultextension='.pdb',
                                       initialdir=os.getcwd(),
                                       filetypes=[("pdb","*.pdb"),("All files","*.*")])
            pdbentry.setvalue(filename)
            return
        def getligand():
            self.ligandfile = tkFileDialog.askopenfilename(defaultextension='.pdb',
                                       initialdir=os.getcwd(),
                                       filetypes=[("mol2","*.mol2"),("All files","*.*")])
            
        def submit():
                  
            #if calcmenu.getcurselection() == 'both':
            #    calcs = ['stability','binding']
            if calcmenu.getcurselection() == 'pka':
                calcs = ['scan']
            else:
                calcs = [calcmenu.getcurselection()]
            mutationlist = mutlist.getvalue().split('\n')
            mutationlist.remove('')
            pdbfile=None; pdb = None
            quality = mutqualentry.getvalue()
            
            if not hasattr(self.DB.meta, 'refprotein') or self.DB.meta.refprotein == None:
                tkMessageBox.showinfo('No ref protein',
                                      'Set a reference (wt) protein first')
                return
            #if self.useref.get() == 1:
            #we use ref pdb by default now
            pdbfile = self.writetempPDB()
            pdbname = self.getrefPDBName()
            
         
            if len(mutationlist) == 0 or mutationlist==[u'']:
                print 'mutation list is empty'
                return
            if hasattr(self.DB.meta,'peatsa_jobs') and nameentry.getvalue() in self.DB.meta.peatsa_jobs:
                print 'job name already used'
                return
            name=nameentry.getvalue()
            expcol = expcolmenu.getcurselection()
            self.submitJob(name=name, pdbname=pdbname,
                           pdb=pdb, pdbfile=pdbfile,
                           ligandfile=self.ligandfile,
                           mutations=mutationlist,
                           calcs=calcs, mutationquality=quality,
                           meta={'expcol':expcol,'pdbname':pdbname})         
            close()
            
        jobdlg = Toplevel()
        jobdlg.geometry('+220+220')
        jobdlg.title('Create Calculation')
        balloon = Pmw.Balloon(jobdlg)
        nameentry = Pmw.EntryField(jobdlg,
                labelpos = 'w',
                label_text = 'Name:',
                validate = validatename,
                value = 'mycalc')
        nameentry.pack(fill=BOTH,expand=1)
        balloon.bind(nameentry, 'Calculation name can be anything, but should be unique')        
        expcols = ['']+self.DB.getSimpleFields()
        expcolmenu = Pmw.OptionMenu(jobdlg,
                labelpos = 'w',
                label_text = 'Exp. col:',
                items = expcols,
                initialitem = '',
                menubutton_width = 8)        
        expcolmenu.pack(fill=BOTH,expand=1)
        balloon.bind(expcolmenu, 'Field with experimental data to compare, optional')          
        calcmenu = Pmw.OptionMenu(jobdlg,
                labelpos = 'w',
                label_text = 'Calculation Type:',
                items = self.calctypes,
                initialitem = 'stability',
                menubutton_width = 8)
        calcmenu.pack(fill=X,expand=1)
        fr=Frame(jobdlg)
        fr.pack(fill=X,expand=1)
        mutqualentry = Pmw.EntryField(jobdlg,
                labelpos = 'w',
                label_text = 'Quality:',
                validate = validatename,
                value = '2.0')
        mutqualentry.pack(fill=BOTH,expand=1)        
        Label(jobdlg,text='Using PDB: '+self.getrefPDBName()).pack(fill=BOTH,expand=1)
        self.ligandfile=None
        mutlist = Pmw.ScrolledText(jobdlg,
                labelpos = 'n',
                label_text='Mutations:',
                usehullsize = 1,
                hull_width = 200,
                hull_height = 250,
                text_wrap='word') 
        mutlist.pack(fill=BOTH,expand=1)
        Button(jobdlg,text='Load Mutations from Project',command=loadmutsfromDB).pack(fill=X,expand=1)        
        Button(jobdlg,text='Load Mutations from File',command=loadmuts).pack(fill=X,expand=1)
        balloon.bind(mutlist, 'Enter one mutation per line in the form\n A:0003:ALA or A3A')
        f=Frame(jobdlg); f.pack(fill=X,expand=1)
        Button(f,text='Submit',command=submit).pack(side=LEFT,fill=X,expand=1,pady=2)
        Button(f,text='Cancel',command=close).pack(fill=X,expand=1,pady=2)        
        jobdlg.grab_set()
        jobdlg.transient(self.parent)
        self.parent.wait_window(jobdlg)
        return
    
    def submitJob(self, name='mycalc', pdbname=None, pdb=None, pdbfile=None, ligandfile=None,
                     mutations=[], calcs=['stability'], mutationquality='2.0', meta={}):
        """Submit job to server"""

        if 'scan' in calcs and pdbname==None:
            print 'You must provide pdb code for pKa calcs'
            return
        if pdb==None and pdbfile==None:
            return
        job = self.jobManager.createJob(pdbId=pdbname, calculations=calcs, 
                                          dataTable='Data', metadata=meta,
                                          optionArgs={'--mutationQuality':mutationquality})
        if pdb != None:
            job.setStructure(pdb)
        else:
            job.setStructureFromFile(pdbfile)
        if 'binding' in calcs:
            job.setLigandFromFile(ligandfile)
        self.mutationList = Core.Data.MutationListFile(filename='tempmutlist', create=True)
        sets=[]
        for code in mutations:
            if code == '': continue
            try:
                sets.append(Core.Data.MutationSet(code))
            except:                
                print 'mutation code %s incorrect' %code
                    
        for s in sets:            
            self.mutationList.addMutant(s, autoUpdate=False, ignoreDuplicates=True)
        self.mutationList.removeDuplicates(autoUpdate=False)
        job.setMutationListFile(self.mutationList)
        job.setState('Ready')
        self.jobManager.logJobStates('jobstates.log')
        #add job to peat database
        self.storeJob(name, job)
        if self.parent != None:
            username = self.parent.username
            self.updateJobs()
        else:
            username = None
        self.DB.commit(note='peatsa job',user=username)
        print 'job submitted successfully'
        return

    def resubmitJob(self):
        """Resend a job based on new mutations in DB that are not in job already"""
        job, name = self.getJob()
        if job == None:
            return
        DB=self.DB       
        self.matrices = job.data.allMatrices()
        for m in matrices:
            matrix=matrices[m]
            if matrix==None: return
            muts = matrix.mutationCodes()
            dbmuts = [DB.get(p).Mutations for p in DB.getRecs()]
            newmuts = list(set(dbmuts) - set(muts))
            print 'the following mutations in the project are not in the job: %s' %newmuts
            
        '''self.submitJob(name=name,
                         pdb=pdb, pdbfile=pdbfile,
                         ligandfile=self.ligandfile,
                         mutations=newmuts,
                         calcs=calcs, meta={'expcol':expcol}) '''
        self.log.yview('moveto', 1)        
        return
        
    def getJob(self, name=None):
        """Get job from name"""
        if name == None:                         
            name = self.jobstable.get_selectedRecordNames()[0]
            if name == None:
                return None, name
        jobid = self.DB.meta.peatsa_jobs[name]
        try:
            job = PEATSA.WebApp.Data.Job(jobid, self.connection)
        except:
            #print 'job not in database'
            return None,name
        return job, name
        
    def removeJob(self):
        """Remove a job from the db"""
        job, name = self.getJob()        
        answer = tkMessageBox.askyesno("Warning",'Remove this job?')
        if answer == False:
            return        
        try:            
            self.jobManager.deleteJob(job)
        except:
            print 'job not in database, removing from peat'
        del self.DB.meta.peatsa_jobs[name]
        self.DB.meta.__p__changed = 1
        self.updateJobs()
        return
        
    def viewDetails(self, name=None):
        job, name = self.getJob()
        if job==None:            
            return
        jobmeta = job.metadata()
        print
        print job.data
        print 'details for job %s' %name
        print 'job status:',job.state()       
        print 'submitted on ',job.date
        if jobmeta.has_key('pdbname'):
            print 'original pdb file:', jobmeta['pdbname']
        print 'mutations:', len(job.mutationListFile().mutantList())
        print '(this job has id %s)' %job.identification        
        if job.error() != None:
            print 'The job had an error..'            
            print job.error()['ErrorDescription']
            print job.error()['DetailedDescription']
        print
        self.log.yview('moveto', 1)
        return

    def addColoredText(self, st, tag, word, fg='black', bg='white'):
        """add a space to the end of the word"""
        word = word + " "
        st.insert('end', word)
        end_index = st.index('end')
        begin_index = "%s-%sc" % (end_index, len(word) + 1)
        st.tag_add(tag, begin_index, end_index)
        st.tag_config(tag, foreground=fg, background=bg)
        return        

    def checkJobsDict(self):
        """Check jobs data structure exists"""
        if not hasattr(self.DB.meta,'peatsa_jobs'):
            from ZODB.PersistentMapping import PersistentMapping
            self.DB.meta.peatsa_jobs = PersistentMapping()
            
    def storeJob(self, name, job):
        """Store job to DB"""
        self.checkJobsDict()
        self.DB.meta.peatsa_jobs[name] = job.identification
        return
            
    def updateJobs(self):
        if not hasattr(self.DB.meta,'peatsa_jobs'):
            return         
        self.updateJobsTable()
        self.wait=self.mainwin.after(60000, self.updateJobs)
        return

    def mergeResults(self, job, colname, tablemodel):
        """Merge given job results to tablemodel"""       
        if job==None:
            return      
        matrices = job.data.allMatrices()
        if not colname:
            return
        nf={'Total':colname}
        for m in matrices:
            matrix = matrices[m]
            if matrix == None: continue
            M = self.mergeMatrix(matrix, tablemodel, fields=['Total'], newfields=nf)            
        return
    
    def mergeCurrent(self):
        """Auto merge selected job results to main table
            called from GUI """
        job, name = self.getJob()
        if job==None:            
            return      
        
        #get field name to use
        colname = tkSimpleDialog.askstring("Column name?",
                                   "Name for column:",
                                    initialvalue=name+'_Predictions',
                                    parent=self.mainwin)
        M = self.parent.tablemodel
        self.mergeResults(job, colname, M)
        self.parent.updateTable()
        
        #also send some meta data to peatsa_meta?
        '''from Correlation import CorrelationAnalyser  
        C = CorrelationAnalyser()    
        cc,rmse = C.getStats(pre,exp)        
        data.append({'name':p,'rmse':rmse,'cc':cc}) '''        
        return
    
    def manageResults(self, name=None):
        """Get the results back - we can send the matrix to the main peat
           table or put results into a labbook sheet.
           Also allow user to merge with an existing table"""
        job, name = self.getJob(name)
        if job.error() != None:
            print 'job had an error, use view details'
        elif job.state() == 'Finished':
            self.showPEATSAResultsDialog(job, name)
        else:
            print 'Job is not finished yet.'
        return            
        
    def editConfigFile(self):
        """Edit config file"""
        from PEATDB.textFrame import textFrame
        tf = textFrame(parent=self.mainwin,
                         title='PEATSA Conf file')
        tf.load_from_file(self.confpath)
        self.parent.wait_window(tf.frame)
        #reconnect
        configuration = Core.Environment.Configuration(filename=self.confpath)
        self.connect(configuration)
        return
    
    def showPEATSAResultsDialog(self, job, name):
        resdlg = Toplevel()
        resdlg.geometry('600x450+300+200')
        resdlg.title('PEATSA results '+name)
        balloon = Pmw.Balloon(resdlg)
        self.currname = name
        body = Frame(resdlg)
        resdlg.initial_focus = body
        body.pack(fill=BOTH,expand=1,padx=5, pady=5) 

        self.matrices = job.data.allMatrices()
        fr=Frame(body)
        fr.grid(row=0,column=0,sticky='news',rowspan=2)
        for m in self.matrices:
            if self.matrices[m] != None:
                self.showMatrix(fr,self.matrices[m], m)
        self.labboklist = self.parent.labbookSheetsSelector(body)  
        self.labboklist.grid(row=0,column=1,sticky='news')                
        bf=Frame(body)
        bf.grid(row=1,column=1,sticky='ew')
        b=Button(bf,text='Merge into main table', command=lambda: self.mergeTable(main=True))
        b.pack(fill=X,expand=1)  
        balloon.bind(b,'Merge results into main DB table')        
        b=Button(bf,text='Merge into Selected', command=self.mergeTable)
        b.pack(fill=X,expand=1)  
        balloon.bind(b,'Merge results into an existing labbook table by matching the mutations')
        b=Button(bf,text='Create new table', command=self.send2Labbook)
        b.pack(fill=X,expand=1)
        balloon.bind(b,'Send results to a new sheet in the main labbook')
        Button(bf,text='Save as CSV', command=self.saveCSV).pack(fill=X,expand=1)      
        body.columnconfigure(0,weight=1)
        body.rowconfigure(0,weight=1)
        return 
        
    def showMatrix(self, frame, matrix, label=''):
        """Show matrix in table"""    
        M = self.matrix2Table(matrix)
        mtable = self.showTable(frame, M, label)
        return mtable
        
    def showTable(self, frame, model, label=''):
        """Show model in table"""
        tf=LabelFrame(frame,text=label)
        tf.pack(fill=BOTH,expand=1)        
        mtable = TableCanvas(tf, model=model, cellwidth=70,                                  
                                 editable=False)
        mtable.createTableFrame()        
        return mtable
        
    def mergeTable(self, main=False):
        """Send a matrix to the peat main table or labbook sheet 
           by merging matching mutations.
           Requires that one field in the table stores compatible 
           mutant format supported by PEATSA"""
        if main == False:
            try:
                name = self.labboklist.getcurselection()[0]
            except:
                print 'no name selected'
                return

        if main == True:
            for m in self.matrices:            
                matrix = self.matrices[m]
                if matrix == None: continue            
                M = self.parent.tablemodel
                M = self.mergeMatrix(matrix, M)
                self.parent.updateTable()
        else:
            M = self.DB.getLabbookSheet(name)
            for m in self.matrices:            
                matrix = self.matrices[m]
                if matrix == None: continue            
                M = self.mergeMatrix(matrix, M)
                if M != None:
                    self.DB.createLabbookSheet(name, M)
                    self.parent.startLabbook('ALL')
        return
        
    def send2Labbook(self):
        """Send matrix to selected labbook"""
        #get name
        cols = ['']+self.DB.getSimpleFields()
        DB=self.DB
        mpDlg = MultipleValDialog(title='Send to Labbook',
                                    initialvalues=(self.currname, cols),
                                    labels=('table name','exp data column'),
                                    types=('string','list'),
                                    parent=self.mainwin)
        if mpDlg.result == False:
            return        
        name = mpDlg.results[0]
        expcol = mpDlg.results[1]
                   
        M = DBActions.sendDB2Labbook(DB,recs=None,cols=['Mutations',expcol],name=name)        
        for m in self.matrices:
            matrix = self.matrices[m]
            if matrix != None:              
                M = self.mergeMatrix(matrix, M)
                self.DB.createLabbookSheet(name, M)
        self.parent.startLabbook('ALL')
        return        
        
    def saveCSV(self):
        """Save matrix to csv"""
        filename=tkFileDialog.asksaveasfilename(defaultextension='.csv',
                                       initialdir=os.getcwd(),
                                       filetypes=[("csv","*.csv"),("All files","*.*")])
        if not filename:
            return
        for m in self.matrices:
            matrix = self.matrices[m] 
            if matrix != None: 
                c=matrix.csvRepresentation()
                f=open(filename,'w')
                f.write(c)
                f.close()
        return
                
    def matrix2Table(self, matrix):
        """Creates a table model from a peatsa matrix"""         
        M = TableModel()
        M.addColumn('Mutations')

        fields = matrix.columnHeaders()
        for f in fields:
            M.addColumn(f)
        i = matrix.indexOfColumnWithHeader('Mutations')
        for row in matrix:
            mutationSet = Core.Data.MutationSet(row[i])
            code = '+'.join(mutationSet.mutationCodes(reduced=True))
            M.addRow(code)
            for f in fields:
                j = matrix.indexOfColumnWithHeader(f)
                if f == 'Mutations':
                    M.data[code]['Mutations'] = code
                else:                   
                    M.data[code][f] = str(row[j])
        return M

    def mergeMatrix(self, matrix, model, fields=None, newfields=None):
        """Merge a peatsa matrix with a table, returns merged tablemodel
        tablemodel: input tablemodel       
        fields: which fields from matrix should be included in merge, default all
        newfields: a dict that can map matrix names to new col names
        """
        M = self.matrix2Table(matrix)
        if fields==None:
            fields = M.columnNames
        key = 'Mutations'        
        if not key in model.columnNames:
            print 'this table has no mutations column, we cannot merge'
            return
        i = matrix.indexOfColumnWithHeader(key)
      
        for row in model.reclist:
            try:
                mset1 = Core.Data.MutationSet(model.data[row][key])
            except:
                continue
            for rec in M.reclist:
                try:
                    mset2 = Core.Data.MutationSet(M.data[rec][key])
                except:
                    continue                 
                if mset1 == mset2:
                    #add this data to table
                    for f in fields:                       
                        if newfields!=None and newfields.has_key(f):
                            col = newfields[f]
                        else:
                            col = f
                        if not M.data[rec].has_key(f): continue
                        model.addColumn(col)
                        try:
                            model.data[row][col] = float(M.data[rec][f])
                        except:
                            model.data[row][col] = M.data[rec][f]
        return model
        
    def showAllResults(self):
        """Show results for single or multiple jobs together"""
        
        names = self.jobstable.get_selectedRecordNames()
        if len(names)==1:
            ax,mh,x,y=self.showResults()

        else:
            tx=[]; ty=[]
            import pylab as plt
            f=plt.figure(figsize=(8,8))
            ax=f.add_subplot(111)
            for n in names:
                a,mh,x,y = self.showResults(n,showtable=False, ax=ax,stats=False)
                tx.extend(x)
                ty.extend(y)                
            ax.legend()
            #add stats for summary
            from Correlation import CorrelationAnalyser
            C = CorrelationAnalyser()           
            C.addStats(ax,tx,ty)
            f.show()
            
        return
        
    def showResults(self, name=None, showtable=True, ax=None, stats=True):
        """Show results with correlation plot from selected job"""
        job, name = self.getJob(name)
        
        if job == None:
            print 'job not in DB'
            return
        if job.state() != 'Finished':
            print 'job not finished'
            return

        self.matrices = job.data.allMatrices()
        #print self.matrices['ModellingResults'].csvRepresentation()
        jobmeta = job.metadata()
        cols = self.DB.getSimpleFields()
        expcol = None
        expdata = None
        #print jobmeta
        if jobmeta.has_key('expcol'):
            expcol = jobmeta['expcol']
        if expcol not in cols and jobmeta.has_key('project'):
            #we may have stored the exp data in another project
            prjdata = jobmeta['project']
            print 'trying to loading exp data from external project(s)'
            from PEATDB.Base import PDatabase
            from PEATTables import PEATTableModel
            
            tmpdb = PDatabase(**prjdata)
            print tmpdb
            S = PEATTableModel(tmpdb)
            expdata = S.simpleCopy(include=['Mutations'])
            print expdata  
            
        #if exp column not known then ask user              
        if expcol == '' or expcol == None:                      
            mpDlg = MultipleValDialog(title='Select Experimental Data',
                                        initialvalues=[cols],
                                        labels=['exp data column:'],
                                        types=['list'],
                                        parent=self.mainwin)
            if mpDlg.result == True:
                expcol = mpDlg.results[0]
            else:
                return

        for m in self.matrices:
            matrix = self.matrices[m]
            if matrix == None or not 'Total' in matrix.columnHeaders():
                continue
            
            ax,mh,x,y = self.plotMerged(matrix, expcol, expdata, m,
                                    showtable, ax, name, stats)
            
            #need to add this for mousehandler to work.. hack       
            '''from Correlation import MouseHandler
            mh = MouseHandler(ax, labels=expcol, key='Mutations')
            mh.connect()'''

        return ax,mh,x,y

    def plotMerged(self, matrix, expcol, expdata=None,
                    title='', showtable=True, ax=None, name=None,
                    stats=True):
        """Merge a set of exp vals with predictions and plot"""
        if expdata==None:
            expdata = self.parent.tablemodel.simpleCopy(include=['Mutations'])
        merged = self.mergeMatrix(matrix, expdata)
        x,y,names,muts = merged.getColumns(['Total',expcol,'name','Mutations'],allowempty=False)
        from Correlation import CorrelationAnalyser
        C = CorrelationAnalyser()
        muts = ['mutation: '+i for i in muts]
        labels = zip(names, muts)
        ax,frame,mh = C.plotCorrelation(x,y,labels,title=title,ylabel=expcol,
                                        ax=ax,plotname=name,stats=stats,err=4)
        x=[round(float(i),2) for i in x]
        y=[round(float(i),2) for i in y]       
        if showtable == True:
            table = self.showTable(frame, merged)
            mh.table = table
            
        return ax,mh,x,y
        
    def test(self):
        job, name = self.getJob('myjob')
        if job.error() != None or job.state() != 'Finished':
            return                                    
        stabmatrix = job.data.stabilityResults
        L = self.DB.getLabbookSheet('myjob')  
        L = self.mergeMatrix(stabmatrix, L, fields=['name'])
        print L.columnNames
        #L1 = self.DB.getLabbookSheet('myjob3')
        #L.merge(L1) 
        return        

    def displayNoDBWarning(self):
        """Warn user that no DB is present"""        
        tkMessageBox.showwarning("Cannot launch plugin",
                         'No Database is currently open. '
                         'You should first open a project.')       
        return
    
    def help(self):
        import webbrowser
        link='http://enzyme.ucd.ie/main/index.php/PEAT_SA'
        webbrowser.open(link,autoraise=1)
        return
        
    def quit(self, evt=None):
        """We MUST stop the jobManager"""
        self.log2Stdout()
        self.jobManager.stopLogging()        
        self.mainwin.destroy()
        print 'closing plugin'
        return  
        
def main():
    import os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="file",
                        help="Open a local db")
    opts, remainder = parser.parse_args()
    #test
    if opts.file != None and os.path.exists(opts.file):
        path=os.path.abspath(opts.file)        
        from PEATDB.Base import PDatabase
        DB = PDatabase(local=path) 
        P = PEATSAPlugin()
        P.main(DB=DB)
        P.test()        
         
if __name__ == '__main__':
    main()

