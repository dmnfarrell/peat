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

import os, sys, copy
import math, random, string, types
import re
import numpy as np
from datetime import datetime
import ConfigParser, csv
from itertools import izip, chain, repeat
import Importer, Custom
from Processing import Processor
import Utilities
from PEATDB.Ekin.Base import EkinProject, EkinDataset
import PEATDB.Ekin.Fitting as Fitting

class Pipeline(object):
    """This class does all the pipeline processing and configuration"""

    __version__ = '1.2.0'
    configsections = ['base','files','fitting','models','variables','functions',
                      'excel','plotting','custom']

    def __init__(self, conffile=None):

        homepath = os.path.join(os.path.expanduser('~'))
        self.defaultpath = os.path.join(homepath, '.pipeline')
        if not os.path.exists(self.defaultpath):
            os.mkdir(self.defaultpath)
        if conffile==None:
            confpath = os.path.join(self.defaultpath, 'default.conf')
            self.createConfig(confpath)
        self.savedir = os.getcwd()
        self.filename = ''
        self.lines = None
        self.queue = {}
        self.results = []
        self.sep = '__'   #symbol for internal separator
        return

    def createConfig(self, conffile='default.conf', **kwargs):
        """Create a basic config file with default options and/or custom values"""

        c = ConfigParser.ConfigParser()
        wdir = os.path.join(self.defaultpath,'workingdir')
        functionsconf = os.path.join(self.defaultpath,'functions.conf')
        defaults = {'base': [('format', 'databyrow'), ('rowstart', 0), ('colstart', 0), ('rowend', 0),
                        ('colend', 0), ('colheaderstart', 0),('rowheaderstart', 0),
                        ('rowheader', ''), ('colheader', ''),
                        ('colheaderlabels', ''), ('rowheaderlabels', ''),
                        ('rowrepeat', 0), ('colrepeat', 0), ('delimeter', ','),
                        ('workingdir', wdir), ('functionsconf',functionsconf),
                        ('ignorecomments', 1),
                        ('decimalsymbol', '.'), ('preprocess',''),
                        ('xformat',''), ('yformat',''), ('groupbyfields', 0)],
                    'files': [('groupbyname', 0), ('parsenamesindex', ''),
                              ('filenameseparator',''), ('parsemethod','words'),
                              ('replicates',0), ('extension','.txt')],
                    'models': [('model1', '')], 'variables': [('variable1', '')],
                    'functions':[('function1','')],
                    'excel': [('sheet', 0), ('numsheets', 1)],
                    'plotting': [('saveplots', 0), ('fontsize',9),
                            ('normalise', 0), ('grayscale', 0), ('alpha',0.8),
                            ('font','sans-serif'), ('markersize',25), ('linewidth',1),
                            ('showerrorbars',0), ('dpi', 100), ('marker','o'),
                            ('markers','-,o'),('legend',0)],
                    'custom': [],
                    'fitting': [('xerror', 0), ('yerror', 0), ('iterations', 50),
                                ('modelsfile','')],
                    }

        cp = Utilities.createConfigParserfromDict(defaults, self.configsections ,**kwargs)
        cp.write(open(conffile,'w'))
        self.parseConfig(conffile)
        return cp

    def parseConfig(self, conffile=None):
        """Parse the config file"""
        f = open(conffile,'r')
        cp = ConfigParser.ConfigParser()

        try:
            cp.read(conffile)
        except Exception,e:
            print 'failed to read config file! check format'
            print 'Error returned:', e
            return
        self.configurationfile = conffile

        #get format so we can create the right Importer
        format = cp.get('base', 'format')
        self.importer = self.getImporter(format, cp)
        if self.importer == None:
            print 'failed to get an importer'
        #both the Importer and Pipeline object get copies of the config options
        #as attributes, convenient but probably bad..
        Utilities.setAttributesfromConfigParser(self, cp)

        self.models = sorted(self.models)
        self.variables = sorted(self.variables)
        print 'parsed config file ok, format is %s' %format
        return

    def writeConfig(self, filename=None):
        """Save a config file from the current object"""
        if filename == None:
            filename = self.configurationfile
        data = self.__dict__
        cp = Utilities.createConfigParserfromDict(data, self.configsections)
        cp.write(open(filename,'w'))
        return

    def getImporter(self, format, cp):
        """Get the required importer object"""
        import inspect
        #try to find a custom importer first
        customclasses = inspect.getmembers(Custom, inspect.isclass)
        defaultclasses = inspect.getmembers(Importer, inspect.isclass)

        for c in customclasses+defaultclasses:
            name, cls = c
            if not hasattr(cls, 'name') or name == 'BaseImporter':
                pass
            elif cls.name == format:
                print 'found importer %s' %name
                return cls(cp)

        print 'no importer found!'
        return None

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

    def loadModels(self):
        #print self.modelsfile
        if self.modelsfile != '':
            try:
                Fitting.loadModelsFile(self.modelsfile)
            except:
                pass
        #print Fitting.currentmodels.keys()
        return

    def openRaw(self, filename=None, callback=None):
        """Open raw file, display preview and get some info about them"""

        if not os.path.exists(filename):
            return None

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

        #override some settings for excel
        self.importer.delimeter = ' '
        return lines

    def doImport(self, lines=None):
        """Import file with current setting and return a dict"""

        if lines == None:
            if self.lines != None:
                lines=self.lines
            else:
                print 'no file loaded'
                return None

        try:
            data = self.importer.doImport(lines)
        except Exception, e:
            print 'Importer returned an error:'
            print e
            #debugging
            import traceback
            tr = sys.exc_info()[2]
            traceback.print_tb(tr)
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

    def prepareData(self):
        """Prepare for data import and processing"""

        if self.workingdir == '':
            self.workingdir = os.path.join(os.getcwd(), 'workingdir')
        if not os.path.exists(self.workingdir):
            os.mkdir(self.workingdir)
        else:
            print 'clearing working directory..'
            Utilities.clearDirectory(self.workingdir)
        return

    def doProcessingStep(self, data, fname):
        """Apply a pre-defined processing step to the data"""

        X = Processor(self.functionsconf)
        names = [i[1] for i in self.functions]
        for n in names:
            if n not in X.predefined:
                print 'function %s not found' %n
                return data
        #user may want to see raw data before processing, so we plot here
        if self.saveplots == 1:
            Er = Utilities.getEkinProject(data)
        data = X.doFunctions(names, data)
        if self.saveplots == 1:
            Em = Utilities.getEkinProject(data)
            Er.addProject(Em,label='proc')
            Er.saveProject(fname+'_processed')
            self.saveEkinPlotstoImages(Er, fname+'_processed',groupbylabel=True)
        return data

    def parseLabels(self):
        """Get labels from filenames"""

        #parseindex can be a list, we currently just take first element
        if type(self.parsenamesindex) is not types.IntType and self.parsenamesindex!='':
            index = self.parsenamesindex.split(',')
            index = int(index[0])
        else:
            index = self.parsenamesindex
        self.namelabels = Utilities.parseFileNames(self.queue, ind=index,
                                                     sep=self.filenameseparator,
                                                     match=self.parsemethod)
        return

    def run(self, callback=None):
        """Do initial import/fitting run with the current config"""

        self.stop=False
        self.loadModels()
        self.prepareData()
        print 'processing files in queue..'

        self.parseLabels()
        imported = {}   #raw data
        results = {}    #fitted data
        #print self.queue

        for key in self.queue:
            filename = self.queue[key]
            lines = self.openRaw(filename)
            if lines == None:
                continue
            data = self.doImport(lines)
            imported[key] = data

        #rebuild dict into a nested structure if it's flat (i.e. from single files)
        '''from Data import NestedData
        D = NestedData(imported)
        D.buildNestedStructure([0,2])
        D.show()
        imported = D.data
        self.namelabels = None'''

        #try to average replicates here before we process
        if self.replicates == 1:
            if self.namelabels != None:
                imported = Utilities.addReplicates(imported, self.namelabels)
            else:
                print 'no replicates detected from labels'

        #re-arrange the imported dict if we want to group our output per field
        if self.groupbyfields == 1:
            imported = Utilities.arrangeDictbySecondaryKey(imported, self.namelabels)

        total = len(imported)
        #print imported
        #print self.namelabels

        c=0.0
        for key in imported:
            if self.stop == True:
                print 'cancelled'
                return
            #set filename
            fname = os.path.basename(key)
            fname = os.path.join(self.workingdir, fname)
            data = imported[key]

            if self.function1 != '':
                data = self.doProcessingStep(data, fname)

            if self.namelabels == None or not self.namelabels.has_key(key):
                namelabel = key
            else:
                namelabel = self.namelabels[key]
            #print namelabel, key
            #print data

            #if we have models to fit this means we might need to propagate fit data
            if self.model1 != '':
                Em = EkinProject()
                #grouping by file labels handled here
                if self.groupbyname == 1:
                    #we don't pass the last model if it has to be
                    #reserved for a final round of fitting from the files dict
                    models = self.models[:-1]
                    variables = self.variables[:-1]
                    E,fits = self.processFits(rawdata=data, Em=Em,
                                               models=models,variables=variables)
                else:
                    E,fits = self.processFits(rawdata=data, Em=Em)
                results[namelabel] = fits
                #print E.datasets, namelabel
            else:
                #if no fitting we just put the data in ekin
                Em = Utilities.getEkinProject(data)
                results[namelabel] = data

            Em.saveProject(fname)
            Em.exportDatasets(fname, append=True)
            if self.model1 != '':
                self.saveFitstoCSV(Em, fname)
            if self.saveplots == 1:
                self.saveEkinPlotstoImages(Em, fname)
            c+=1.0
            if callback != None:
                callback(c/total*100)

        #if grouped by file names then we process that here from results
        if self.groupbyname == 1:
            results = Utilities.extractSecondaryKeysFromDict(results)
            Em = EkinProject()
            #print results
            E,fits = self.processFits(rawdata=results, Em=Em)
            fname = os.path.join(self.workingdir, 'final')
            Em.saveProject(os.path.join(self.workingdir, fname))
            Em.exportDatasets(os.path.join(self.workingdir, fname))
            if self.model1 != '':
                self.saveFitstoCSV(Em, fname)
            self.saveEkinPlotstoImages(Em, fname)
        print 'processing done'
        print 'results saved to %s' %self.workingdir
        self.results = results
        return results

    def processFits(self, rawdata, models=None, variables=None, ind=None,
                    parentkey='', Em=None):
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
        #print nesting,currmodel,currvariable,parentkey

        if nesting == 0:
            #bottom level of nesting, we just fit
            xerror = float(self.xerror); yerror = float(self.yerror)
            E = Utilities.getEkinProject(rawdata, xerror=xerror, yerror=yerror)
            print rawdata.keys(), parentkey, currmodel
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
                #now we pass each child node recursively
                E,fit = self.processFits(rawdata[l], ind=ind-1, parentkey=lbl, Em=Em)
                fitdata[l] = fit
            E = Utilities.getEkinProject(fitdata)
            if parentkey == '': parentkey = 'final'
            E,fit = self.getFits(E, currmodel, currvariable, str(parentkey))
            Em.addProject(E,label=parentkey)
            return E,fit

    def getFits(self, E, model, varname='a', filename=None):
        """Fit an Ekin project
           model: model to fit
           varname: variable to extract"""

        fits = []
        xerrors = []
        yerrors = []
        E.fitDatasets('ALL', models=[model], noiter=self.iterations,
                      conv=1e-9, grad=1e-8, silent=True)
        labels = E.datasets
        if self.xerror != 0 or self.yerror != 0:
            print 'getting exp uncert'
            self.expUncert(E, xerr=float(self.xerror))
            for d in labels:
                yerrors.append(E.getMeta(d,'exp_errors')[varname][1])

        for d in labels:
            fits.append(E.getMetaData(d)[varname])
        if self.saveplots == 1 and filename != None and filename != '':
            print 'plotting %s' %filename
            self.saveEkinPlotstoImages(E, filename)
        return E,(labels,fits,xerrors,yerrors)

    def findVariableName(self, model, i=0):
        """Find a parameter name for the given model at index i,
           used so we don't cause an error when doing processFits when no
           variable is given in the conf file"""

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
        """Find if bottom level of nesting or not"""
        for d in data.keys():
            if type(data[d]) is types.DictType:
                return 1
            else:
                return 0

    def addFolder(self, path):
        """Add a folder to the queue"""

        D = {}
        ext = self.extension
        #print os.path.basename(path), os.path.normpath(path)
        for root, dirs, files in os.walk(path):
            for d in dirs:
                if d.startswith('.'):
                    dirs.remove(d)
            for f in files:
                fname = os.path.join(root,f)
                if ext in os.path.splitext(fname)[1] or ext == '':
                    #print root,fname
                    label = os.path.relpath(fname,  os.path.normpath(path))
                    label = label.replace(os.path.sep,'_')
                    #print label
                    D[label] = fname
        self.queue = D
        return D

    def addtoQueue(self, files):
        """Add files"""

        if type(files) is types.StringType or type(files) is types.UnicodeType:
            files = [files]
        for f in files[:]:
            if f not in self.queue.values():
                self.queue[f] = f
        return

    def getPlotOptions(self):
        """Get options dict for passing to plotdatasets"""
        opts = dict(self.plotting)
        #print self.plotting
        for k in opts:
            if k=='markers':
                opts[k]=list(opts[k].split(','))
            else:
                try:
                    opts[k] = eval(opts[k])
                except:
                    pass
        return opts

    def saveEkinPlotstoImages(self, E, filename, groupbylabel=False):
        """Save the ekin plots to png images
           groupbylabel allows commonly labelled datasets
           to be overlayed - useful for viewing processing results"""
        title = os.path.basename(filename)
        filename = os.path.join(self.workingdir, filename)
        filename = filename + '.png'
        #limit to 30 plots per image
        dsets = E.datasets
        if len(dsets) > 30:
            dsets = E.datasets[:30]
        plotopts = self.getPlotOptions()
        #the following functionality should be integrated into ekin plotter
        if groupbylabel == True:
            import itertools
            import pylab as plt
            from mpl_toolkits.axes_grid1 import AxesGrid
            groups = []
            #use iterator to get groups
            iterator = itertools.groupby(sorted(dsets), key=lambda x: x.split('_')[0])
            for i,j in iterator:
                groups.append(list(j))
            #print groups
            fig = plt.figure()
            cols, dim = E.getGeometry(groups, 0)
            c=1
            for g in groups:
                ax=fig.add_subplot(cols,dim,c)
                E.plotDatasets(g,figure=fig,axis=ax,plotoption=3,
                                varyshapes=1,**plotopts)
                c+=1
            #fig.tight_layout(pad=0.8, w_pad=0.8, h_pad=1.0)
            fig.savefig(filename,dpi=plotopts['dpi'])
        else:
            E.plotDatasets(dsets, filename=filename, plotoption=2,
                           size=(10,8), title=title, **plotopts)
        return

    def saveFitstoCSV(self, E, filename):
        """Save results to a csv file"""
        title = os.path.basename(filename)
        filename = os.path.join(self.workingdir, filename)
        filename = filename + '_fits.csv'
        from PEATDB.Ekin.Tables import EkinProjModel
        EM = EkinProjModel(E)
        EM.exportCSV(filename)
        return

