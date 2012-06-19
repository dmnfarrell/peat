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

    __version__ = '1.0.0'
    configsections = ['base','files','fitting','models','variables','functions',
                      'excel','plotting','custom']

    def __init__(self, conffile=None):

        if conffile==None:
            confpath = os.path.join(os.path.expanduser('~'),'default.conf')
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
        wdir = os.path.join(os.path.expanduser('~'),'workingdir')
        functionsconf = os.path.join(os.getcwd(),'functions.conf')
        defaults = {'base': [('format', 'databyrow'), ('rowstart', 0), ('colstart', 0), ('rowend', 0),
                        ('colend', 0), ('colheaderstart', 0),('rowheaderstart', 0),
                        ('rowheader', ''), ('colheader', ''),
                        ('rowrepeat', 0), ('colrepeat', 0), ('delimeter', ','),
                        ('workingdir', wdir), ('functionsconf',functionsconf),
                        ('ignorecomments', 1),
                        ('decimalsymbol', '.'), ('preprocess',''),
                        ('xformat',''),('yformat','')],
                    'files': [('groupbyname', 0), ('parsenamesindex', 0),
                                ('replicates',0)],
                    'models': [('model1', '')], 'variables': [('variable1', '')],
                    'functions':[('function1','')],
                    'excel': [('sheet', 0), ('numsheets', 1)],
                    'plotting': [('saveplots', 0), ('fontsize',9),                                    
                            ('normalise', 0), ('grayscale', 0), ('alpha',0.8),
                            ('font','sans-serif'), ('markersize',25), ('linewidth',1),
                            ('showerrorbars',0), ('dpi', 100), ('marker','o')],
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
        print self.modelsfile
        if self.modelsfile != '':
            try:
                Fitting.loadModelsFile(self.modelsfile)
            except:
                pass
        print Fitting.currentmodels.keys()
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
            Er = self.getEkinProject(data)            
        data = X.doFunctions(names, data)
        if self.saveplots == 1:
            Em = self.getEkinProject(data)      
            Er.addProject(Em,label='proc')
            Er.saveProject(fname+'_processed')
            self.saveEkinPlotstoImages(Er, fname+'_processed',groupbylabel=True)
        return data

    def run(self, callback=None):
        """Do initial import/fitting run with the current config"""

        self.stop=False
        self.loadModels()
        self.prepareData()
        print 'processing files in queue..'

        if self.groupbyname == 1:
            labels = self.parseFileNames(self.queue, ind=self.parsenamesindex)
        else:
            labels = self.queue

        self.labels = labels
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

        total = len(imported)
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
            Em.exportDatasets(fname)
            if self.model1 != '':
                self.saveFitstoCSV(Em, fname)

            if self.saveplots == 1:
                self.saveEkinPlotstoImages(Em, fname)
            c+=1.0
            if callback != None:
                callback(c/total*100)

        #if grouped by file names then we process that here from results
        if self.groupbyname == 1:
            #print results.keys()
            results = self.extractSecondaryKeysFromDict(results)
            Em = EkinProject()
            E,fits = self.processFits(rawdata=results, Em=Em)
            fname = os.path.join(self.workingdir, 'final')
            Em.saveProject(os.path.join(self.workingdir, fname))
            Em.exportDatasets(os.path.join(self.workingdir, fname))
            if self.model1 != '':
                self.saveFitstoCSV(Em, fname)
            #if self.saveplots == 1:
            self.saveEkinPlotstoImages(Em, fname)
        print 'processing done'
        print 'results saved to %s' %self.workingdir
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
                arrs.append(np.array(D[n]))
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
            print ind, f, labels[f]
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
            if opts[k] == '1':
                opts[k] = True
            elif opts[k] == '0':
                opts[k] = False
            if k=='alpha':
                opts[k]=float(opts[k])
            if k in ['dpi','markersize']:
                opts[k]=int(opts[k])
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
                E.plotDatasets(g,figure=fig,axis=ax,plotoption=3,markers=['-','o'],
                                varyshapes=1,**plotopts)
                c+=1
            fig.tight_layout(pad=0.8, w_pad=0.8, h_pad=1.0)
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

    @classmethod
    def getEkinProject(self, data, xerror=None, yerror=None, sep='__'):
        """Get an ekin project from a dict of the form
             {label:([x],[y]),..} or
             {label:([x],[y],[xerr],[yerr]),..}"""

        E = EkinProject(mode='General')
        for d in data.keys():
            if type(data[d]) is types.DictType:
                for lbl in data[d]:
                    name = str(d)+sep+str(lbl)
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


