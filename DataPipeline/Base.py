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

import os, sys, math, random, numpy, string, types
from datetime import datetime
import ConfigParser, csv

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
        
