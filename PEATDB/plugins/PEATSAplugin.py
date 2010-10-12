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
import tkFileDialog
import Pmw
import PEATSA.WebApp.Data
import PEATSA.WebApp.UtilityFunctions
import PEATSA.Core as Core
from PEATDB.Dialogs import MultipleValDialog
from PEATDB.Actions import DBActions
from PEATDB.TableModels import TableModel
from PEATDB.Tables import TableCanvas
import tkMessageBox

class PEATSAPlugin(Plugin):
    """Template GUI plugin for PEAT App"""
    capabilities = ['gui']
    requires = ['PEATSA']
    menuentry = 'PEATSA Plugin'
    gui_methods = {'createJobDialog':'Create New Job',
                   'fetchJob':'Fetch Job from Server',                  
                   'editConfigFile' : 'Configure Server',
                   'help':'Help',
                   'quit':'Close'}
    buttonorder = ['createJobDialog','fetchJob','editConfigFile','help','quit']
    about = 'This plugin allows you to call PEATSA'    

    calctypes = ['stability','binding','both']
    
    def main(self, parent=None, DB=None):       
        if parent==None:
            if DB!=None:
                self.DB = DB
            else:
                return 
        else:
            self.parent = parent
            self.DB = parent.DB
            self.parentframe = None            
            self._doFrame()
            self.setupConnection()
            print 'Updating jobs table..'
            self.updateJobs()            
        return

    def setupConnection(self):
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
            if parent != None:               
                tkMessageBox.showwarning("Connection Error",
                         'No PEATSA server configured, press configure server'
                         ' to set a server, username and password.')
        self.connect(configuration)
        return
    
    def _doFrame(self):
        self.mainwin = self.parent.createChildFrame(width=450,title='PEATSA Plugin')        
        methods = self._getmethods()
        methods = [m for m in methods if m[0] in self.gui_methods.keys()]        
        l=Label(self.mainwin, text='Run PEAT-SA Calculations')
        l.pack(side=TOP,fill=BOTH)
        self.tf=LabelFrame(self.mainwin,text='Jobs in this project')
        self.tf.pack(side=TOP,fill=BOTH,expand=1)       
        self.manageJobsButtons(self.mainwin)    
        self._createButtons(methods)
        self.log = self.createLogWin(self.mainwin)
        self.log.pack(side=TOP,fill=BOTH,expand=1)         
        self.stdout2Log()
        return

    def _createButtons(self, methods):
        """Dynamically create buttons for supplied methods, which is a tuple
            of (method name, label)"""     
        for m in methods:
            b=Button(self.mainwin,text=self.gui_methods[m[0]],command=m[1])
            b.pack(side=BOTTOM,fill=BOTH)
        return

    def updateJobsTable(self):
        """Show table for current jobs list"""        
        self.checkJobsDict()
        jobdict = self.DB.meta.peatsa_jobs
        M = TableModel()
        #open job log from file
        jl = pickle.load(open('jobstates.log','r'))
        for j in jobdict: 
            jobid =jobdict[j]           
            try:                
                M.addRecord(j,state=jl[jobid]['State'],date=jl[jobid]['Date'])
            except:
                M.addRecord(j,state='Not in DB')
        self.jobstable = TableCanvas(self.tf, model=M, height=100, editable=False)
        self.jobstable.createTableFrame()       
        self.log.yview('moveto', 1)        
        return

    def manageJobsButtons(self, parent):
        fr1 = Frame(parent)
        Button(fr1,text='View Results',command=self.showResults,bg='#ccFFFF').pack(side=TOP,fill=BOTH,expand=1)
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
            import tkFileDialog
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
                  
            if calcmenu.getcurselection() == 'both':
                calcs = ['stability','binding']
            else:
                calcs = [calcmenu.getcurselection()]
            mutationlist = mutlist.getvalue().split('\n')
            mutationlist.remove('')
            pdbfile=None; pdb = None
            
            if not hasattr(self.DB.meta, 'refprotein') or self.DB.meta.refprotein == None:
                tkMessageBox.showinfo('No ref protein',
                                      'Set a reference (wt) protein first')
                return     
            #if self.useref.get() == 1:
            #we use ref pdb by default now
            pdbfile = self.writetempPDB()
               
            '''else:
                p = pdbentry.getvalue()
                if os.path.exists(p):
                    pdbfile=p
                elif len(p) == 4:                    
                    pdb = DBActions.fetchPDB(p)       
                else:
                    print 'no valid pdb given'
                    return  '''        
         
            if len(mutationlist) == 0 or mutationlist==[u'']:
                print 'mutation list is empty'
                return
            if hasattr(self.DB.meta,'peatsa_jobs') and nameentry.getvalue() in self.DB.meta.peatsa_jobs:
                print 'job name already used'
                return
            name=nameentry.getvalue()
            expcol = expcolmenu.getcurselection()
            self.submitJob(name=name,
                           pdb=pdb, pdbfile=pdbfile,
                           ligandfile=self.ligandfile,
                           mutations=mutationlist,
                           calcs=calcs, meta={'expcol':expcol})         
            close()
            
        jobdlg = Toplevel()
        jobdlg.geometry('+220+220')
        jobdlg.title('Create new job')
        balloon = Pmw.Balloon(jobdlg)
        nameentry = Pmw.EntryField(jobdlg,
                labelpos = 'w',
                label_text = 'Job name:',
                validate = validatename,
                value = 'myjob')
        nameentry.pack(fill=BOTH,expand=1)
        balloon.bind(nameentry, 'Job name can be anything, but should be unique')        
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
                label_text = 'Calculation:',
                items = self.calctypes,
                initialitem = 'stability',
                menubutton_width = 8)
        calcmenu.pack(fill=X,expand=1)
        fr=Frame(jobdlg)
        fr.pack(fill=X,expand=1)
        '''self.useref = IntVar()
        Label(fr,text='Use reference protein').pack(side=LEFT)
        Checkbutton(fr,variable=self.useref).pack(side=LEFT)        
        pdbentry = Pmw.EntryField(jobdlg,
                labelpos = 'w',
                label_text = 'PDB Structure:')
        pdbentry.pack(fill=X,expand=1)
        balloon.bind(pdbentry, 'Enter the PDB ID or load a file')
        Button(jobdlg,text='load PDB from file',command=getstruct).pack(fill=X,expand=1)
        Button(jobdlg,text='load Ligand file',command=getligand).pack(fill=X,expand=1)'''
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
    
    def submitJob(self, name='myjob', pdb=None, pdbfile=None, ligandfile=None,
                     mutations=[], calcs=['stability'], meta={}):    
        """Submit job to server"""

        if pdb==None and pdbfile==None:
            return
        job = self.jobManager.createJob(name, calculations=calcs, 
                                          dataTable='Data', metadata=meta)
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
        self.DB.commit(note='peatsa job',user=self.parent.username)  
        self.updateJobs()
        print 'job submitted successfully'
        return

    def resubmitJob(self):
        """Resend a job based on new mutations in DB that are not in job already"""
        job, name = self.getJob()
        if job == None:
            return
        DB=self.DB
        dataset = job.data
        matrices = {'binding':dataset.bindingResults,
                    'stability':dataset.stabilityResults}
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
        print
        print 'details for job %s' %name
        print 'job status:',job.state()
        #self.addColoredText(self.log,'blu', 'job status: '+job.state() , fg='blue')
        print 'submitted on ',job.date
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
    
    def help(self):
        import webbrowser
        link='http://enzyme.ucd.ie/main/index.php/PEAT_SA'
        webbrowser.open(link,autoraise=1)
        return
    
    def quit(self):      
        self.mainwin.destroy()
        self.jobManager.stopLogging()
        self.log2Stdout()
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
        dataset = job.data
        self.matrices = {'binding':dataset.bindingResults,
                         'stability':dataset.stabilityResults}
        fr=Frame(body)
        fr.grid(row=0,column=0,sticky='news',rowspan=2)
        for m in self.matrices:
            if self.matrices[m] != None:
                self.showMatrix(fr,self.matrices[m], m)
        self.labboklist = self.parent.labbookSheetsSelector(body)
        #self.labboklist.insert(0,'main db')
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
        for m in self.matrices:
            matrix = self.matrices[m]
            if matrix == None: continue
            if main == True:
                M = self.parent.tablemodel
                M = self.mergeMatrix(matrix, M)
                self.parent.updateTable()
            else:
                M = self.DB.getLabbookSheet(name)
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

    def mergeMatrix(self, matrix, model, fields=None):
        """Merge a peatsa matrix with a table, returns merged tablemodel
        tablemodel: input tablemodel
        key: use given key to match on, usually mutations code
        fields: which fields from model should be included in merge, default is all
        """
        if fields==None:
            fields = model.columnNames
        key = 'Mutations'
        M = self.matrix2Table(matrix)      
        if not key in model.columnNames:
            print 'this table has no mutations column, we cannot merge'
            return
        i = matrix.indexOfColumnWithHeader(key)
      
        for row in M.reclist:           
            mset1 = Core.Data.MutationSet(M.data[row][key])
            for rec in model.reclist:
                try:
                    mset2 = Core.Data.MutationSet(model.data[rec][key])
                except:
                    continue                 
                if mset1 == mset2:
                    #add this data to table
                    for f in fields:
                        if not model.data[rec].has_key(f): continue
                        #print f, model.data[rec][f]
                        M.addColumn(f)                      
                        M.data[row][f] = model.data[rec][f]
        return M        
        
    def showResults(self):
        """Show results with correlation plot from selected job"""
        job, name = self.getJob()
        if job == None:
            print 'job not in DB'
            return
        if job.state() != 'Finished':
            print 'job not finished'
            return
        dataset = job.data
        self.matrices = {'binding':dataset.bindingResults,
                         'stability':dataset.stabilityResults}
        jobmeta = job.metadata()
        if jobmeta.has_key('expcol'):
            expcol = jobmeta['expcol']
        else:
            expcol=''
        if expcol == '' or expcol==None:
            #if exp column not known then ask user
            cols = self.DB.getSimpleFields()            
            mpDlg = MultipleValDialog(title='Select Experimental Data',
                                        initialvalues=[cols],
                                        labels=['exp data column:'],
                                        types=['list'],
                                        parent=self.mainwin)
            if mpDlg.result == True:
                expcol = mpDlg.results[0]
            else:
                return
        
        #from PEATDB.plugins.Correlation import CorrelationAnalyser  
        from Correlation import CorrelationAnalyser  
        C = CorrelationAnalyser()        
        
        for m in self.matrices:
            matrix = self.matrices[m]          
            if matrix == None: continue
            M = self.parent.tablemodel           
            M = self.mergeMatrix(matrix, M, fields=['name',expcol])           
            x,y,names,muts = M.getColumns(['Total',expcol,'name','Mutations'],allowempty=False)          
            labels = zip(names, muts)
            ax,frame,mh = C.plotCorrelation(x,y,labels,title=m,ylabel=expcol)
            table = self.showTable(frame, M)
            mh.table = table
        return
    
    def test(self):
        job, name = self.getJob('myjob')
        if job.error() != None or job.state() != 'Finished':
            return                                    
        stabmatrix = job.data.stabilityResults
        L = self.DB.getLabbookSheet('myjob')
        #n,m,x,y = L.getColumns(['name','Mutations','Total','deltatm'],allowempty=False)
        #print x,y,t,d        
        L = self.mergeMatrix(stabmatrix, L, fields=['name'])
        print L.columnNames
        #L1 = self.DB.getLabbookSheet('myjob3')
        #L.merge(L1)
 
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

