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
from Importer import *
import Utilities
from PEATDB.Ekin.Base import EkinProject, EkinDataset

class Pipeline(object):
    """This class does all the pipeline processing and configuration"""

    sections = ['base','fitting','excel']

    def __init__(self, conffile=None):

        if conffile==None:
            self.createConfig('default.conf')
        self.savedir = os.getcwd()
        self.filename = ''
        self.lines = None
        self.queue = {}
        self.results = []
        self.sep = '__'   #symbol for internal separator
        return

    def createConfig(self, filename, **kwargs):
        """Create a basic config file with default options and/or custom values"""

        c = ConfigParser.ConfigParser()
        wdir = os.path.join(os.getcwd(),'workingdir')

        defaults = {'base': [('format', 'databyrow'), ('rowstart', 0), ('colstart', 0), ('rowend', 0),
                        ('colend', 0), ('rowheader', ''), ('colheader', ''),
                        ('rowrepeat', 0), ('colrepeat', 0), ('delimeter', ','),
                        ('workingdir', wdir),  ('ignorecomments', 1),
                        ('checkunicode', 0), ('decimalsymbol', '.')],
                    'files': [('groupbyname', 0), ('parsenamesindex', 0),
                                ('replicates',0)],
                    'models': [('model1', '')], 'variables': [('variable1', '')],
                    'excel': [('sheet', 0), ('numsheets', 1)],
                    'plotting': [('saveplots', 0), ('normaliseplots', 0), ('grayscale', 0),
                        ('dpi', 100)],
                    'custom': [],
                    'fitting': [('xerror', 0), ('yerror', 0), ('iterations', 30), ('modelsfile','')],
                    }
        order = ['base','files','fitting','models','variables',
                 'excel','plotting','custom']

        for s in order:
            c.add_section(s)
            for i in defaults[s]:
                name,val = i
                c.set(s, name, val)

        #use kwargs to create specific settings in the appropriate section
        for s in c.sections():
            opts = c.options(s)
            for k in kwargs:
                if k in opts:
                    c.set(s, k, kwargs[k])
        #handle model and variable sections which can have zero or multiple
        #options
        for k in sorted(kwargs):
            if k.startswith('model'):
                c.set('models', k, kwargs[k])
            elif k.startswith('variable'):
                c.set('variables', k, kwargs[k])

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
        #as attributes, convenient but probably bad..
        Utilities.setAttributesfromConfigParser(self, cp)
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
                              datainrows=0,modelsfile='test.dict')
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
        self.addtoQueue(filename)
        self.filename = filename
        self.lines = lines
        return lines

    def closeFile(self):
        """Close current file"""
        self.filename = ''
        self.lines = None
        return

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
                print 'no file loaded'
                return None

        #self.parseConfig(self.conffile)
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

    def preProcess(self):
        """Prepare for data import and processing"""
        if not os.path.exists(self.workingdir):
            os.mkdir(self.workingdir)
        else:
            print 'clearing working directory..'
            Utilities.clearDirectory(self.workingdir)
        return

    def run(self, callback=None):
        """Do initial import/fitting run with the current config"""

        self.preProcess()
        print 'processing files in queue..'

        if self.groupbyname == 1:
            labels = self.parseFileNames(self.queue, ind=self.parsenamesindex)
        else:
            labels = self.queue

        self.labels = labels
        total = len(self.queue)
        c=0.0
        imported = {}   #raw data
        results = {}    #fitted data

        #print 'queue'
        #print self.queue
        #print filelabels

        for key in self.queue:
            filename = self.queue[key]
            lines = self.openRaw(filename)
            data = self.doImport(lines)
            imported[key] = data

        #try to average replicates here before we process
        if self.replicates == 1:
            imported = self.handleReplicates(imported)

        for key in imported:
            #set filename
            fname = os.path.basename(key)
            fname = os.path.join(self.workingdir, fname)
            data = imported[key]
            if labels.has_key(key):
                label = labels[key]
            else:
                label = key
            #print key,fname,label

            #if we have models to fit this means we might need to propagate fit data
            if self.model1 != '':
                Em = EkinProject()
                if self.groupbyname == 1:
                    #we don't pass the last model if it has to be
                    #reserved for a final round of fitting from the files dict
                    models = self.models[:-1]
                    variables = self.variables[:-1]
                    E,fits = self.processFits(rawdata=data, Em=Em,
                                               models=models,variables=variables)
                else:
                    E,fits = self.processFits(rawdata=data, Em=Em)
                results[label] = fits
            else:
                 #if no fitting we just put the data in ekin
                Em = self.getEkinProject(data)
            Em.saveProject(fname)

            if self.saveplots == 1:
                self.saveEkinPlotstoImages(Em, fname)
            c+=1.0
            if callback != None:
                callback(c/total*100)

        #if groupby files then we process that here from results
        if self.groupbyname == 1:
            results = self.extractSecondaryKeysFromDict(results)
            #print results
            Em = EkinProject()
            E,fits = self.processFits(rawdata=results, Em=Em)
            fname = os.path.join(self.workingdir, 'final')

            Em.saveProject(os.path.join(self.workingdir, fname))
            if self.saveplots == 1:
                self.saveEkinPlotstoImages(Em, fname)
        print 'processing done'
        return

    def processFits(self, rawdata, models=None, variables=None, ind=None, parentkey='', Em=None):
        """Process the a dict of possibly nested dicts
            ind: the index indicating the level of recursion, used to find the right model
                 and variable
            parentkey: the label from the parent dict that will be matched to the fits
            returns: final set of fits"""

        nesting = self.getDictNesting(rawdata)
        if models == None:
            models = self.models
            variables = self.variables

        if len(models) == 0:
            print 'no models found for fitting!'
            return None, None
        if ind == None:
            ind = len(models)-1

        def getmodelinfo():
            model = models[ind][1]
            try:
                var = variables[ind][1]
            except:
                var = ''
            if var == '':
                var = self.findVariableName(model)
            return model, var

        currmodel,currvariable = getmodelinfo()
        #print models, variables
        #print nesting,currmodel,currvariable

        if nesting == 0:
            #final level of nesting, we just fit
            xerror = float(self.xerror); yerror = float(self.yerror)
            E = self.getEkinProject(rawdata, xerror=xerror, yerror=yerror)
            E,fit = self.getFits(E, currmodel, currvariable, str(parentkey))
            Em.addProject(E, label=parentkey)
            return E,fit
        else:
            #if there is nesting we pass the subdicts recursively and get their fits
            fitdata = {}
            for l in rawdata.keys():
                if parentkey != '':
                    lbl = str(l)+'_'+str(parentkey)
                else:
                    lbl = l
                #now we pass each child node to the same function
                E,fit = self.processFits(rawdata[l], ind=ind-1, parentkey=lbl, Em=Em)
                fitdata[l] = fit
            E = self.getEkinProject(fitdata)
            if parentkey == '': parentkey = 'final'
            E,fit = self.getFits(E, currmodel, currvariable, str(parentkey))
            Em.addProject(E,label=parentkey)
            return E,fit

    def handleReplicates(self, data):
        """If the configuration specifies replicates than we average the
           corresponding data points in the raw dict """

        print 'processing replicates..'
        import operator
        from itertools import groupby
        newdata = {}
        labels = self.labels
        sorteditems = sorted(labels.iteritems(), key=operator.itemgetter(1))

        for key, group in groupby(sorteditems, lambda x: x[1]):
            c = 0
            subdata = []
            for g in group:
                f = g[0]    #key in corresponding data dict
                c+=1
                subdata.append(data[f])
                print f
            newdata[key] = self.averageDicts(subdata)
            print '%s replicates for label %s' %(c, key)
        #print newdata
        return newdata

    def averageDicts(self, dictslist):
        """Average dicts of the form
           {label1: [[x1],[y1]],..],label2:[[x2],[y2]]...}"""

        newdata = {}
        names = dictslist[0].keys()
        for n in names:
            arrs = []
            for D in dictslist:
                arrs.append(numpy.array(D[n]))
            #print sum(arrs)/len(arrs)
            newdata[n] = list(sum(arrs)/len(arrs))
        return newdata

    def getFits(self, E, model, varname='a', filename=None):
        """Fit an Ekin project
           model: model to fit
           varname: variable to extract"""

        fits = []
        xerrors = []
        yerrors = []
        E.fitDatasets('ALL', models=[model], noiter=self.iterations,
                      conv=1e-6, grad=1e-6, silent=True)
        labels = E.datasets
        if self.xerror != 0 or self.yerror != 0:
            print 'getting exp uncert'
            self.expUncert(E, xerr=float(self.xerror))
            for d in labels:
                yerrors.append(E.getMeta(d,'exp_errors')[varname][1])

        for d in labels:
            fits.append(E.getMetaData(d)[varname])
        if self.saveplots == 1 and filename != None and filename != '':
            self.saveEkinPlotstoImages(E, filename)
        return E,(labels,fits,xerrors,yerrors)

    def findVariableName(self, model, i=0):
        """Find a parameter name for the given model at index i,
           used so we don't cause an error when doing processFits when no
           variable is given in the conf file"""

        import PEATDB.Ekin.Fitting as Fitting
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

    def getDictNesting(self, data):
        """Get level of nesting"""
        for d in data.keys():
            if type(data[d]) is types.DictType:
                return 1
            else:
                return 0

    def parseFileNames(self, filenames, ind=0):
        """Parse file names to extract a numerical value
           ind: extract the ith instance of a number in the filename"""
        labels = {}
        for f in filenames:
            bname = os.path.basename(f)
            l = re.findall("([0-9.]*[0-9]+)", bname)[ind]
            labels[f] = l
            #print f, labels[f]
        return labels

    def parseString(self, filename):
        """Parse string to extract a pattern"""

        return

    def addFolder(self, path, ext='.txt'):
        """Add a folder to the queue"""

        D = {}
        #print os.path.basename(path), os.path.normpath(path)
        for root, dirs, files in os.walk(path):
            for d in dirs:
                if d.startswith('.'):
                    dirs.remove(d)
            for f in files:
                fname = os.path.join(root,f)
                #absname = os.path.abspath(fname)
                if ext in os.path.splitext(fname)[1]:
                    #print root,fname
                    label = os.path.relpath(fname,  os.path.normpath(path))
                    label = label.replace(os.path.sep,'_')
                    #print label
                    D[label] = fname
        #self.addtoQueue(D.values())
        self.queue = D
        return D

    def addtoQueue(self, files):
        """Add files"""
        if type(files) is types.StringType:
            files = [files]
        for f in files[:]:
            if f not in self.queue.values():
                self.queue[f] = f
        return

    def saveEkinPlotstoImages(self, E, filename):
        """Save the ekin plots to png images"""
        title = os.path.basename(filename)
        filename = os.path.join(self.workingdir, filename)
        filename = filename + '.png'
        #limit to 25 plots per image
        d = E.datasets
        if len(d) > 25: d = E.datasets[:25]
        E.plotDatasets(d, filename=filename, plotoption=2,
                       dpi=self.dpi, size=(10,8), showerrorbars=True,
                       title=title, normalise=self.normaliseplots,
                       grayscale=self.grayscale)
        return

    @classmethod
    def getEkinProject(self, data, xerror=None, yerror=None):
        """Get an ekin project from a dict of the form
             {label:([x],[y]),..} or
             {label:([x],[y],[xerr],[yerr]),..}"""

        E = EkinProject(mode='General')
        for d in data.keys():
            if type(data[d]) is types.DictType:
                for lbl in data[d]:
                    name = d+sep+lbl
                    xy = data[d][lbl]
                    ek=EkinDataset(xy=xy)
                    E.insertDataset(ek, name)
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
    def extractSecondaryKeysFromDict(self, data):
        """Re-arrange a dict of the form {key1:([names],[values]),..
           into the form {name1:([keys],[values]),..}"""

        newdata = {}
        kys = data.keys()
        datasets = data[kys[0]][0]
        for d in datasets:
            newdata[d] = []
        for l in data:
            names, vals, xerr, yerr = data[l]
            #print zip(names, vals)
            for d,v in zip(names, vals):
                newdata[d].append((l,v))
        for d in newdata:
            newdata[d] = zip(*newdata[d])
        return newdata

class FolderStructure(object):
    """Class to hold structure rawdata that makes handling it easier"""

    def __init__(self, path):
        self.path = path
        self.rootdir = os.path.basename(path)
        self.info = Utilities.getdirectoryStructure(path)
        print self.rootdir
        print '------------------'
        print self.info
        self.traverse(self.info)
        return

    def traverse(self, data):
        for i in data:
            if type(data[i]) is types.DictType:
                self.traverse(data[i])
            else:
                print i, data[i]

    def report(self):
        """report on structure found in path"""
        print 'rootdir is %s' %self.rootdir
        folders = len(self.info)
        print '%s files found in %s folders' %(self.p)
        return

