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

import os, sys
import math, random, numpy, string, types
import re
from datetime import datetime
import ConfigParser, csv
from itertools import izip, chain, repeat
from PEATDB.Ekin.Base import EkinProject, EkinDataset

class Pipeline(object):
    """This class does all the pipeline processing and configuration"""

    sections = ['base','fitting','excel']

    def __init__(self, conffile=None):

        if conffile==None:
            self.createConfig('default.conf')
        self.savedir = os.getcwd()
        self.filename = ''
        self.lines  =None
        self.queue = []
        self.results = []
        self.sep = '__'   #symbol for internal separator
        return

    def createConfig(self, filename, **kwargs):
        """Create a basic config file with default options and/or custom values"""

        c = ConfigParser.ConfigParser()
        s = 'base'
        c.add_section(s)
        c.set(s, 'format', 'databyrow')
        c.set(s, 'delimeter', ',')
        c.set(s, 'decimalsymbol', '.')
        c.set(s, 'checkunicode', 0)
        wdir = os.path.join(os.getcwd(),'workingdir')
        c.set(s, 'workingdir', wdir)
        c.set(s, 'rowstart', 0)
        c.set(s, 'colstart', 0)
        c.set(s, 'rowend', 0)
        c.set(s, 'colend', 0)
        #c.set(s, 'alternatecols', 0)
        #c.set(s, 'alternaterows', 0)
        c.set(s, 'rowrepeat', 0)
        c.set(s, 'colrepeat', 0)
        #c.set(s, 'rowdataformat', 'number')
        #c.set(s, 'coldataformat', 'number')
        c.set(s, 'timeformat', "%Y-%m-%d %H:%M:%S")
        c.set(s, 'ignorecomments', 1)
        c.set(s, 'rowheader', '')
        c.set(s, 'colheader', '')
        c.set(s, 'groupbyfile', 0)
        #fitting settings
        f = 'fitting'
        c.add_section(f)      
        c.set(f, 'xerror', 0)
        c.set(f, 'yerror', 0)
        c.set(f, 'iterations', 30)        
        #models
        m = 'models'
        c.add_section(m)
        c.set(m, 'model1', '')
        c.set(m, 'model2', '')
        c.set(m, 'model3', '')
        #variables
        v = 'variables'
        c.add_section(v)
        c.set(v, 'variable1', '')
        c.set(v, 'variable2', '')
        c.set(v, 'variable3', '')
        #excel settings
        e='excel'
        c.add_section(e)
        c.set(e, 'sheet', 0)
        c.set(e, 'numsheets', 1)

        m='custom'
        c.add_section(m)

        #use kwargs to create specific settings in the appropriate section
       
        for s in c.sections():
            opts = c.options(s)
            #print s, c.options(s)
            for k in kwargs:
                if k in opts:
                    #print s,k
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

        #get format so we can create the right Importer
        format = cp.get('base', 'format')
        self.importer = self.getImporter(format, cp)
        if self.importer == None:
            print 'failed to get an importer'
        #both the Importer and Pipeline object get copies of the config options
        #as attributes
        self.setAttributesfromConfigParser(self, cp)
        self.models = self.getListFromConfigItems(self.models)
        self.variables = self.getListFromConfigItems(self.variables)
        print 'parsed config file ok, format is %s' %format
        return
        
    def getImporter(self, format, cp):
        """Get the required importer object"""
        import Custom
        import inspect
        #try to find a custom importer first
        customclasses = inspect.getmembers(Custom, inspect.isclass)
        for f in Custom.importers:
            if f == format:
                classname = Custom.importers[f]
                importerclass = getattr(Custom, classname)
                print 'found custom importer %s' %classname
                return importerclass(cp)

        if format == 'databyrow':
            importer = DatabyRowImporter(cp)
        elif format == 'databycolumn':
            importer = DatabyColImporter(cp)
        elif format == 'paireddatabyrow':
            importer = PairedDatabyRowImporter(cp)
        elif format == 'paireddatabycolumn':
            importer = PairedDatabyColImporter(cp)
        elif format == 'groupeddatabyrow':
            importer = GroupedDatabyRowImporter(cp)
        elif format == 'groupeddatabycolumn':
            importer = GroupedDatabyColImporter(cp)
        else:
            importer = None
            print 'no importer found!'
        return importer

    def loadPreset(self, preset=None):
        """Load preset config for specific type data"""
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
        if not filename in self.queue:
            self.queue.append(filename)
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
            print 'Importer returned an error:'
            print e
            data = {}
        self.checkImportedData(data)
        return data

    def checkImportedData(self, data):
        """Report info on imported dict"""
        
        l=0
        for d in data.keys():
            if type(data[d]) is types.DictType:
                l+=len(data[d].keys())
            else: l+=1
        print 'importer returned %s datasets' %l
        return
        
    def run(self, callback=None):
        """Do initial import/fitting run with the current config"""

        if not os.path.exists(self.workingdir):
            os.mkdir(self.workingdir)

        print 'processing files in queue.'
        if self.groupbyfile == True:
            filelabels = self.parseFileNames(self.queue)
        else:
            filelabels = None

        total = len(self.queue)
        c=0.0
        #rawdata = self.rawdata = {}
        results = {}
        for filename in self.queue:
            if filelabels != None:
                key = filelabels[filename]
            else:
                key = filename
            lines = self.openRaw(filename)            
            data = self.doImport(lines)
            print self.model1
            #if we have models to fit this means we might need to propagate fit data          
            if self.model1 != '' and self.groupbyfile == 0:
                #pass an ekin project to store all the fits for later viewing
                Em = EkinProject()  
                E,fits = self.processFits(rawdata=data, Em=Em)
                print 'final fits', fits
                results[key] = fits
                #save fits/plots to ekin
                fname = os.path.basename(filename)
                fname = os.path.join(self.workingdir, fname+'.ekinprj')
                Em.saveProject(fname)
                print Em, fname

            c+=1.0
            if callback != None:
                callback(c/total*100)
                
        #if groupbyfiles then we process that here from results, 
        #but labels come from filenames....        
        #e.g.
        #if self.groupbyfile == 1:
            #fits = self.getFits(results, nextmodel, filename)
        
        print 'initial processing done'

        return

    def getDictNesting(self, data):
        """Get level of nesting"""
        for d in data.keys():
            if type(data[d]) is types.DictType:
                return 1
            else:
                return 0
                
    def processFits(self, rawdata=None, ind=None, parentkey='', Em=None):
        """Process the a dict of possibly nested dicts
            ind: the index indicating the level of recursion, used to find the right model
                 and variable
            parentkey: the label from the parent dict that will be matched to the fits
            returns: final set of fits"""
            
        models = self.models
        variables = self.variables
        
        if len(models) == 0:
            print 'no models found for fitting!'
            return None, None
        if ind == None:
            ind = len(models)-1            
        currmodel = models[ind]
        if len(variables) == 0:            
            currvariable = self.findVariableName(currmodel)
            print 'no variable given to extract, using %s' %currvariable
        else:
            currvariable = variables[ind]

        nesting = self.getDictNesting(rawdata) 
        
        #print ind, models, variables
        #print nesting, models[ind]

        if nesting == 0:
            #final level of nesting, we just fit
            xerror = float(self.xerror); yerror = float(self.yerror)
            E = self.getEkinProject(rawdata, xerror=xerror, yerror=yerror)
            E,fit = self.getFits(E, currmodel, currvariable)
            Em.addProject(E, label=parentkey)
            return E,fit
        else:
            #if there is nesting we pass the subdicts recursively and get their fits 
            fitdata = {}
            for l in rawdata.keys():
                print l
                if parentkey!='':
                    lbl = str(l)+'_'+str(parentkey)
                else:
                    lbl = l
                #now we pass each child node to the same function    
                E,fit = self.processFits(rawdata[l], ind=ind-1, parentkey=lbl, Em=Em)
                fitdata[l] = fit
  
            E = self.getEkinProject(fitdata)
            E,fit = self.getFits(E, currmodel, currvariable)
            #print fit
            Em.addProject(E,label=parentkey)
            return E,fit
    
    def getFits(self, E, model, varname='a'):
        """Fit an Ekin project   
           model: model to fit
           varname: variable to extract"""
           
        fits = []
        xerrors = []
        yerrors = []
        E.fitDatasets('ALL', models=[model], noiter=20,
                      conv=1e-6, grad=1e-6, silent=True)
        labels = E.datasets
        if self.xerror != 0 or self.yerror != 0:
            print 'getting exp uncert'
            self.expUncert(E, xerr=float(self.xerror))
            for d in labels:
                yerrors.append(E.getMeta(d,'exp_errors')[varname][1])
        for d in labels:
            fits.append(E.getMetaData(d)[varname])            
        return E,(labels,fits,xerrors,yerrors)
        
    def findVariableName(self, model, i=0):
        """Find a parameter name for the given model at index i,
           used so we don't cause an error when doing processFits when no
           variable is given in the conf file"""
            
        from PEATDB.Ekin.Fitting import Fitting    
        X = Fitting.getFitter(model)
        if X==None:
            print 'no fitter found for this model name'
            return 'a'
        varname = X.getVarNames()[i]
        return varname
        
    def expUncert(self, E, xerr=0, yerr=0):
        """Get fitted variable errors using the iterative method"""
    
        for d in E.datasets:
            ferrs = E.estimateExpUncertainty(d, runs=10, xuncert=xerr, yuncert=yerr)
            E.addMeta(d, 'exp_errors', ferrs)
        E.saveProject()
        return

    def parseFileNames(self, filenames):
        """Parse file names to extract a numerical value"""
        labels = {}
        for f in filenames:
            bname = os.path.basename(f)
            l = re.findall("([0-9.]*[0-9]+)", bname)[0]
            labels[f] = l
        return labels
        
    def addFolder(self, path, ext='.txt'):
        """Add a folder to the queue"""

        for root, dirs, files in os.walk(path):
            for d in dirs:
                if d.startswith('.'):
                    dirs.remove(d)
            #print root, dirs, files
            for f in files:
                fname = os.path.join(root,f)
                if os.path.splitext(fname)[1] == ext:
                    self.addtoQueue([fname])
        #print self.queue
        #self.parseFileNames(self.queue)
        return

    def addtoQueue(self, files):
        """Add files"""
        for f in files[:]:
            if f not in self.queue:
                self.queue.append(f)
        return

    @classmethod
    def getEkinProject(self, data, xerror=None, yerror=None):
        """Get an ekin project from a dict of the form
             {label:([x],[y]),label2:([x],[y])} or
             {label:([x],[y],[xerr],[yerr]),label2:([x],[y],[xerr],[yerr])}"""

        E = EkinProject(mode='General')
        for d in data.keys():
            if type(data[d]) is types.DictType:
                pass
            else:
                #print data[d]
                if len(data[d]) == 4:
                    x,y,xerrs,yerrs = data[d]
                else:
                    x,y = data[d]
                    xerrs = []; yerrs=[]
                    if xerror!=None:
                        xerrs=[xerror for i in x]
                    if yerror!=None:
                        yerrs=[yerror for i in y]                   
                ek = EkinDataset(xy=[x,y], xerrs=xerrs, yerrs=yerrs)
                E.insertDataset(ek, d)                
                #print ek.errors
        return E

    @classmethod
    def setAttributesfromConfigParser(self, obj, cp):
        """A helper method that makes the options in a ConfigParser object
           attributes of obj"""
       
        for s in cp.sections():
            obj.__dict__[s] = cp.items(s)          
            #print cp.items(s)
            for f in cp.items(s):
                #print f[0], f[1]
                try: val=int(f[1])
                except: val=f[1]
                obj.__dict__[f[0]] = val
                
    @classmethod
    def getListFromConfigItems(self, items):
        """Get a list from a set of ConfigParser key-value pairs"""
        lst = [i[1] for i in items if i[1] != '']
        return lst
        
class BaseImporter(object):
    """Base Importer class, sub-class this to define methods specific to each kind of
       import format. At minimum we override the doImport method to get specific
       functionality"""

    def __init__(self, cp):
        """Arguments:
            cp - a ConfigParser object that has been loaded in the parent app"""
        Pipeline.setAttributesfromConfigParser(self, cp)
        if self.delimeter=='': self.delimeter=' '
        elif self.delimeter=='tab': self.delimeter='\t'
        return

    def guessRowStart(self, lines):
        """If rowstart is not specified in config, it might be non-zero"""
        s = 0
        for line in lines:
            if not self.checkEmptyRow(line):
                s+=1
            else:
                self.rowstart = s
                return

    def checkEmptyRow(self, line):
        if line.startswith('#') or line=='' or line.startswith('\r'):
            return False
        return True

    def getRow(self, lines, row, grouped=False):
        """Return values in a row"""
        #if we have repeating cols we return multiple lists
        if self.ignorecomments==True and lines[row].startswith('#'):
            return None
        if not self.checkEmptyRow(lines[row]):
            return None
        vals = string.strip(lines[row]).split(self.delimeter)
        vals = vals[self.colstart:]
        if grouped == False:
            return vals
        else:
            if self.colrepeat == 0:
                return [vals]
            elif self.colrepeat > 1:
                return self.groupList(self.colrepeat, vals)
            else:
                return None

    def getColumn(self, lines, col, grouped=False):
        """Return values in a column"""
        vals = []
        for row in range(self.rowstart, self.rowend):
            if self.ignorecomments==True and lines[row].startswith('#'):
                continue
            rowdata = string.strip(lines[row]).split(self.delimeter)
            vals.append(rowdata[col])
        if grouped == False:
            return vals
        else:
            if self.rowrepeat == 0:
                return [vals]
            elif self.rowrepeat > 1:
                return self.groupList(self.rowrepeat, vals)

    def getColumnHeader(self, lines, grouped=False):
        """Column headers are taken from colstart row"""
        if self.rowheader == '':
            row = self.rowstart
        else:
            row = self.rowheader
        return self.getRow(lines, row, grouped)

    def getRowHeader(self, lines, grouped=False):
        """Return values in header row"""
        if self.colheader == '':
            col = self.colstart
        else:
            col = self.colheader
        return self.getColumn(lines, col, grouped)

    def groupList(self, n, l, padvalue=None):
        """group a list into chunks of n size"""
        return [l[i:i+n] for i in range(0, len(l), n)]

    def getXYValues(self, xd, yd, start=0):
        """Return a lists of floats from lists of vals whilst doing
           various checks to remove errors etc."""

        x=[]; y=[]
        #print xd, yd
        for i in range(start,len(xd)):
            if i>=len(yd): break
            xval = self.checkValue(xd[i])
            yval = self.checkValue(yd[i])
            if xval==None or yval==None:
                continue
            x.append(xval)
            y.append(yval)
        if len(x)<=1:
            return None
        return x,y

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
    """This importer handles data formatted in columns with common x values
       in a specified column, it also handles multiple sets of data grouped
       in evenly spaced distances down each column"""
    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):

        data = {}
        if self.rowend == 0:
            self.rowend=len(lines)
        xdata = self.getRowHeader(lines,grouped=True)
        header = self.getColumnHeader(lines)
        if self.colend == 0:
            self.colend = len(header)
        if xdata == None:
            return

        for col in range(self.colstart+1, self.colend):
            coldata = self.getColumn(lines, col, grouped=True)
            #print xdata, coldata
            if coldata == None: continue
            for xd,yd in zip(xdata,coldata):
                if len(xd)<=1 or len(yd)<=1: continue
                name = yd[0]
                x,y = self.getXYValues(xd[1:],yd[1:])
                data[name] = [x,y]
        return data

class DatabyRowImporter(BaseImporter):
    """This importer handles data formatted in rows with common x values
       along the top (or specified in a specific row), it also handles
       multiple sets of data grouped in evenly spaced distances along the row"""
    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):

        data = {}
        self.guessRowStart(lines)
        xdata = self.getColumnHeader(lines,grouped=True)
        if xdata == None:
            return
        if self.rowend == 0:
            self.rowend=len(lines)

        for row in range(self.rowstart+1, self.rowend):
            if row>=len(lines):
                break
            rowdata = self.getRow(lines,row,grouped=True)
            #print xdata, rowdata
            if rowdata == None: continue
            for xd,yd in zip(xdata,rowdata):
                #print xd, yd
                if len(xd)<=1 or len(yd)<=1: continue
                name = yd[0]
                x,y = self.getXYValues(xd[1:],yd[1:])
                data[name] = [x,y]
        return data

class PairedDatabyColImporter(BaseImporter):
    """This importer handles data formatted in rows with paired x-y values,
       there are therefore no common x-values in the header """
    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):

        data = {}
        if self.rowend == 0:
            self.rowend=len(lines)
        header = self.getColumnHeader(lines)
        if self.colend == 0:
            self.colend = len(header)

        for col in range(self.colstart+1, self.colend):
            coldata = self.getColumn(lines, col, grouped=True)
            for xyd in coldata:
                name = xyd[0]
                x = xyd[1:len(xyd):2]
                y = xyd[2:len(xyd):2]
                data[name] = [x,y]

        return data

class PairedDatabyRowImporter(BaseImporter):
    """This importer handles data formatted in rows with paired x-y values,
       there are therefore no common x-values in the header """
    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):

        data = {}
        if self.rowend == 0:
            self.rowend=len(lines)

        for row in range(self.rowstart+1, self.rowend):
            if row>=len(lines):
                break
            rowdata = self.getRow(lines,row, grouped=True)
            for xyd in rowdata:
                name = xyd[0]
                x = xyd[1:len(xyd):2]
                y = xyd[2:len(xyd):2]
                data[name] = [x,y]

        return data

class GroupedDatabyRowImporter(BaseImporter):
    """This importer handles data formatted in rows with multiple independent x values in
       each column, each dataset is then repeated in groups every x rows, specified in
       the rowrepeat option. The importer therefore returns dictionary with multiple sets of
       x-y values for each label/dataset"""
    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):

        data = {}
        #assumes the column header has labels for each set of xy vals
        labels = self.getColumnHeader(lines)
        print labels
        if self.rowend == 0:
            self.rowend=len(lines)
        if self.colend == 0:
            self.colend = len(labels)

        grouplen = (self.rowend - self.rowstart) / self.rowrepeat
        step = self.rowrepeat

        for d in range(1,grouplen):
            for row in range(self.rowstart, self.rowend, step):
                if row>=len(lines):
                    break
                xdata = self.getRow(lines, row)
                rowdata = self.getRow(lines, row+d)
                name = rowdata[0]
                if not data.has_key(name):
                    data[name] = {}

                for v in zip(labels,xdata,rowdata)[1:]:
                    label = v[0]
                    x = self.checkValue(v[1])
                    y = self.checkValue(v[2])
                    if x==None or y==None:
                        continue
                    if not data[name].has_key(label):
                        data[name][label]=[]
                    l = data[name][label]
                    l.append((x,y))

        #reformat paired vals into x and y lists
        for d in data:
            for lbl in data[d]:
                data[d][lbl] = zip(*data[d][lbl])

        return data

class GroupedDatabyColImporter(BaseImporter):
    """This importer handles data formatted in cols with multiple independent x values in
       each row, each dataset is then repeated in groups every x columns, specified in
       the colrepeat option. The importer therefore returns dictionary with multiple sets of
       x-y values for each label/dataset"""
    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):

        data = {}
        #assumes the row header has labels for each set of xy vals

        if self.rowend == 0:
            self.rowend=len(lines)
        labels = self.getRowHeader(lines)
        header  = self.getColumnHeader(lines)
        if self.colend == 0:
            self.colend = len(header)
        grouplen = len(self.getColumnHeader(lines, grouped=True)[0])
        step = self.colrepeat

        for d in range(1,grouplen):
            for col in range(self.colstart, self.colend, step):

                xdata = self.getColumn(lines, col)
                coldata = self.getColumn(lines, col+d)
                name = coldata[0]
                if not data.has_key(name):
                    data[name] = {}
                #print name, xdata,coldata
                for v in zip(labels,xdata,coldata)[1:]:
                    label = v[0]
                    x = self.checkValue(v[1])
                    y = self.checkValue(v[2])
                    if x==None or y==None:
                        continue
                    if not data[name].has_key(label):
                        data[name][label]=[]
                    l = data[name][label]
                    l.append((x,y))

        #reformat paired vals into x and y lists
        #by convention we put secondary labels into nested dicts
        for d in data:
            for lbl in data[d]:
                data[d][lbl] = zip(*data[d][lbl])
        return data
