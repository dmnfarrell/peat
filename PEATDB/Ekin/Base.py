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

import math, sys, os, types
import IO
from Convert import EkinConvert
from PEATDB.Ekin.Dataset import EkinDataset
from List_Utils import *
from Fitting import Fitting
from Meta import MetaData
import numpy
import Utils
import matplotlib
try:
    matplotlib.use('Agg')
except:
    pass
try:
    from matplotlib.figure import Figure
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
except:
    pass

from Pylab import Options


class EkinProject(object):
    """Ekin base class - capturing, storing and analysing biophysical data.
       This class provides the base functionality that can be called outside
       a GUI.
       Usage: E = EkinProject(data=ekindata)
    """

    ekinfields = ['__datatabs_fits__', '__fit_matches__', '__Exp_Meta_Dat__',
                  '__distance_matches__', '__datatab_structmapping__', '__meta_data__',
                  '__plotopts__', '__currentdataset__', '__displaymultiple__']

    modes =['General',
                'Simple enzyme kinetic',
                'Enzyme pH-activity',
                'pH-stability',
                'NMR titration',
                'Protein Stability',
                'Amyloid formation']
    mode_definition = {'General':
                        ['Michaelis-Menten',#'Michaelis-Menten-competitive-inhibition',
                         'Linear','Power','Gaussian',
                         '1 pKa 2 Chemical shifts','3 pKas, 2 Chemical shifts',
                         '2 pKas, 3 Chemical shifts','3 pKas, 4 Chemical shifts',
                         'Bell-shaped pH-act profile (3 pKas)',
                         'Bell-shaped pH-act profile (2 pKas)',
                         'Sigmoid', 'Modified Hill', 'schellman',
                         'Chemical Denaturation','diffDenaturation','Unfolding',
                         'Residual Activity','Arrhenius',
                         'DSC2state','DSCindependent','DSC2stateIrreversible','DSC2stateIrreversibleII',
                         'Amyloid Fibre Formation'],
                     'Simple enzyme kinetic':['Michaelis-Menten','Michaelis-Menten-competitive-inhibition',
                                              'Competitive inhibition','Non-competitive inhibition'],
                     'Enzyme pH-activity':['Bell-shaped pH-act profile (3 pKas)','Bell-shaped pH-act profile (2 pKas)'],
                     'pH-stability':['Bell-shaped pH-act profile (3 pKas)','Bell-shaped pH-act profile (2 pKas)'],
                     'NMR titration':[ 'Linear', '1 pKa 2 Chemical shifts','3 pKas, 2 Chemical shifts','2 pKas, 3 Chemical shifts',
                                      '3 pKas, 4 Chemical shifts'],#,'4 pKas, 5 Chemical shifts'],
                     'Protein Stability':['Sigmoid'],
                     'Amyloid formation':['Linear','Amyloid Fibre Formation']}


    def __init__(self, parent=None, ekinproj=None, data=None, mode=None,
                  protein=None, field=None):
        """Initialize the class"""

        #try to load pylab if present
        self.protein = protein
        self.parent = parent
        self.currentmode = mode
        #set a default model for fitting new datasets
        if self.currentmode != None:
            self.defaultmodel = self.mode_definition[self.currentmode][0]

        if self.parent and protein:
            self.allowreturntoDB = 1
            self.protein_name=self.parent.data['DBinstance'].DB[protein]['Name']
            #check for other open instances for the cell if called in PEAT
            if self.parent.ekin_instances.has_key(protein+field):
                if self.parent.ekin_instances[protein+field] > 1:
                    print 'EKIN ALREADY OPENED FOR THIS CELL'
                    self.allowreturntoDB = 0
        else:
            self.protein_name=protein
        self.field=field
        self.filename = None
        #create IO vars
        self.savedir = os.getcwd()

        #if ekinprojectfile provided we load it
        if ekinproj != None:
            self.openProject(ekinproj)
        else:
            #just load the data
            self.loadData(data)
        self.current = None    #track 'current' dataset
        self.length = 0        
        return

    def loadData(self, data=None):
        """Load the data into self.data and update selector etc"""
        import copy
        # Did we get data or should we just add an empty data window?
        if data == None or type(data) is types.StringType:
            self.data={}
            self.datasets=[]
            self.__datatabs_fits__= {}
            M = MetaData()
            self.__meta_data__ = M.data
        else:
            self.data = copy.deepcopy(data)
            self.datasets=[]
            # Set mode if present
            if self.data.has_key('mode'):
                self.currentmode = self.data['mode']
                del self.data['mode']
            #we now load the special dicts directly, so same names
            for name in self.ekinfields:
                  if data.has_key(name):
                      self.__dict__[name]=copy.deepcopy(self.data[name])
                      del self.data[name]
            #populate list of datasets
            for dataset in self.data.keys():
                if dataset in self.ekinfields:
                    continue
                self.datasets.append(dataset)
        self.length = len(self.datasets)
        self.checkMeta()
        return

    def checkMeta(self):
        """Load old style meta data if present and then remove"""
        if not hasattr(self,'__meta_data__'):
            M=MetaData(project=self)
            if hasattr(self, '__Exp_Meta_Dat__'):
                print 'old meta data - converting'
                M.convertfromOld()
                del self.__Exp_Meta_Dat__
            self.setMetaData(M.data)
        return

    def addDataset(self, label=None, update=1):
        """Add a new empty dataset, can be combined with insert? """

        if label==None:
            return
        if label in self.datasets:
            print 'error, dataset present'
            return

        self.datasets.append(label)
        self.length = len(self.datasets)
        #here we need to add to self.data
        ek=self.data[label]=EkinDataset()
        for i in range(8):
            ek.add()            
        self.__datatabs_fits__[label]={}
        return

    def insertDataset(self, newdata=None, newname='data', xydata=None, fit=None,
                              replace=False, update=1):
        """Insert pre-existing data as a dataset"""
        #print newdata
        if newdata == None:            
            ek = EkinDataset(xy=xydata)
        else:
            #check for old or bad data
            ek = self.checkDataset(newdata)            
            if ek == False:
                return
                
        if self.data.has_key(newname) and replace==False:
            newname = newname + 'a'        
        self.data[newname] = ek
        self.datasets.append(newname)
        self.length = len(self.datasets)
        if fit==None:
            fit = ek.getFit()
        self.setFitData(newname,fit)
        return

    def readDataset(self, filename):
        """Get dataset from file"""
        newdata = {}
        if os.path.isfile(filename):
            fd=open(filename)
            import pickle
            newdata=pickle.load(fd)
            fd.close()
            # Add the data tab with a filename
            showname=os.path.split(filename)[1]
            showname=showname.replace('_','-')
        return newdata

    def deleteDataset(self, name=None):
        """Delete a dataset"""
        if name == None:
            print 'error, no dataset given'
            return
        if not name in self.datasets:
            print 'error, no such dataset'
            return
        del self.data[name]
        self.datasets.remove(name)
        self.length = len(self.datasets)
        for f in self.ekinfields:
            if not self.__dict__.has_key(f) or type(self.__dict__[f]) is not types.DictType:
                continue
            if self.__dict__[f].has_key(name):
                del self.__dict__[f][name]
        return

    def renameDataset(self, currname, newname):
        """Rename a dataset"""
        if newname == None or newname == currname:            
            return
        self.copyDataset(currname, newname)
        self.deleteDataset(currname)
        return

    def copyDataset(self, name, newname=None):
        """Copy a dataset"""

        if newname == None:
            newname = name + '_copy'
        elif newname in self.datasets:
            print 'error, new name exists'
            return
        else:
            import copy
            self.data[newname] = copy.deepcopy(self.data[name])
            for f in self.ekinfields:
                if not self.__dict__.has_key(f) or type(self.__dict__[f]) is not types.DictType:
                    continue
                if self.__dict__[f].has_key(name):
                    self.__dict__[f][newname] = copy.deepcopy(self.__dict__[f][name])
            self.datasets.append(newname)
            self.length = len(self.datasets)
        return

    def getDataset(self, dataset):
        """Get a dataset"""
        try:
            return self.data[dataset]
        except:
            return None

    def deleteAll(self):
        """Delete All"""
        self.data={}
        self.datasets=[]
        self.__datatabs_fits__={}
        self.length = 0
        return

    def saveDataset(self, name, filename=None):
        """Save single dataset to a file"""
        import tkFileDialog, os
        if filename == None:
            filename = name + '.ekindat'

        ek = self.data[name]
        fitdata = self.getFitData(name)
        ek.addFit(fitdata)
        fd=open(filename,'w')
        import pickle
        pickle.dump(ek,fd)
        fd.close()
        return

    def checkDict(self, dataset, name):
        """Check for arbitrary fields stored in the old format dataset
           we move these to metadata"""
        if not type(dataset) is types.DictType:
            return None        
        for k in dataset.keys():
            if not k in [0,1,2]:            
                self.addMeta(name, k, dataset[k])
        #print self.getMetaData(name)
        return        
        
    def checkDataset(self, dataset):
        """Check for old format or bad dataset"""           
        if type(dataset) is types.DictType:
            #old format
            x,y,a,xerr,yerr = EkinConvert.ekin2xy(dataset, getall=1, geterrors=True)           
            return EkinDataset(x=x,y=y,active=a,xerrs=xerr,yerrs=yerr)           
        elif not type(dataset) is EkinDataset:
            return False
        else:
            return dataset
                
    def checkDatasets(self):
        """Check corrupt datasets and remove"""
        #print 'checking for old or bad datasets'        
        for d in self.datasets[:]:
            ek = self.checkDataset(self.data[d])            
            if ek == False:
                print 'removed', d
                self.deleteDataset(d)
            else:
                other = self.checkDict(self.data[d], d)
                self.data[d] = ek                
        return

    def importCSV(self, filename=None, text=None, format=1, sep=',',
                  mode=None, truncate='end'):
        """
        Import csv file(s) to ekin data.
        This method should work without user intervention as much as possible.
        There are three options:
            1. Single file, curves per row, with ph values at the top
            2. Single file, curves per column, with ph values in
            the first column and curve names in first row
            3. Each curve in its own file, with 2 cols
        """
        import csv

        if filename != None:
            csvfile = open(filename)
            reader = csv.reader(csvfile, delimiter=sep)

        elif text != None:
            reader = csv.reader(text.split(os.linesep))

        x = reader.next()
        pts = len(x)
        print 'found csv data with %s points in header row' %pts
        datasets = {}
        for row in reader:
            if len(row)==0:
                continue
            name=row[0]; del row[0]
            xdata=x
            if len(row) < pts:
                print 'dataset is short %s points' %(pts-len(row))
                if truncate == 'end':
                    xdata = x[0:len(row)]
                else:
                    st=pts-len(row)
                    xdata = x[st:len(x)]
            elif len(row) > pts:
                row = row[0:pts]
                print 'dataset has more than %s points, truncating' %pts

            #datasets[name] = EkinConvert.xy2ekin([xdata,row])
	    datasets[name] = EkinDataset(xy=[xdata,row])
        for d in datasets.keys():
            self.insertDataset(datasets[d], d)
        print 'imported %s datasets' % len(datasets)
        print '-----------------------'
        return

    def exportCSV(self, datasets='ALL', filename=None, sep=','):
        """Export dataset(s) to csv"""
        import csv
        if filename != None:
            csvfile = open(filename, 'w')
        writer = csv.writer(csvfile, delimiter=sep)
        datasets = self.checkDatasetParam(datasets)

        for d in self.datasets:
            writer.writerow([d])            
	    ek=self.getDataset(d)
	    x,y = ek.getxy()
            for i in range(len(x)):
                writer.writerow([x[i],y[i]])
        return

    def openProject(self, filename=None):
        """Open an ekin project file"""
        import os
        if filename==None:
            return
        if os.path.splitext(filename)[1] == '':
            filename = filename + '.ekinprj'        
        if not os.path.exists(filename):
            #print 'no such file'
            return None
        data={}
        if filename != None:
            if os.path.isfile(filename):
                fd=open(filename)
                import pickle
                data=pickle.load(fd)
                fd.close()

        self.loadData(data)
        self.filename = filename       
        return

    def saveProject(self, filename=None):
        """Save data to an ekin project file"""
        if filename == None:
            if not hasattr(self, 'filename') or self.filename==None:
                print 'no filename'
                return
            else:
                filename=self.filename
        self.filename = filename
        if os.path.splitext(filename)[1] == '':
            filename = filename + '.ekinprj'
        data = self.prepare_data()
        fd=open(filename,'w')
        import pickle
        pickle.dump(data, fd)
        fd.close()
        return

    def prepare_data(self):
        """Prepare data for saving or return to another app ie. PEAT"""
        import copy
        data = copy.deepcopy(self.data)
        data = self.data
        #print self.__meta_data__
        for name in self.ekinfields:
            if hasattr(self, name):
                data[name] = copy.deepcopy(self.__dict__[name])
        data['mode'] = self.currentmode

        return data

    def printData(self):
        """Print out a summary of all the datasets"""
        for name in self.datasets:
            fitdata = self.__datatabs_fits__[name]
            if fitdata.has_key('model'):
                fit=fitdata['model']
            else:
                fit='Not fitted'
            axes = self.getAxesLabels(name)
            print '%s: %s pts, fit: %s, axes: %s' %(name,len(self.data[name][0].keys()),fit, axes)

    def prettyPrint(self):
        """prints the entire dict"""
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        x=pp.pformat(self.data)
        print x
        return

    def getFitData(self, dataset):
        """Get the fit data for the dataset"""
        if self.__datatabs_fits__.has_key(dataset):
            return self.__datatabs_fits__[dataset]
        else:
            return None

    def setFitData(self, dataset, fitdata=None):
        """Set the fit data for the dataset"""
        if fitdata != None:
            self.__datatabs_fits__[dataset] = fitdata
        else:
            self.__datatabs_fits__[dataset] = {}
        return fitdata

    def addMeta(self, dataset, fieldname, data):
        """Set a field for this dataset in the meta data dict"""
        if not hasattr(self, '__meta_data__'):
            self.__meta_data__ = {}
            self.__meta_data__[dataset] = {}
        if not  self.__meta_data__.has_key(dataset):
            self.__meta_data__[dataset] = {}
        self.__meta_data__[dataset][fieldname] = data
        return

    def getMetaData(self, dataset):
        """Return the fit data in a formatted dict with the variable
           names as keys, more intelligible than the simple fit data.
           We also include the experimental meta data here, if present
           and all fields which are special keys e.g. __enzyme_conc__"""
        mdata={}
        fitdata = self.__datatabs_fits__[dataset]
        if fitdata != None and len(fitdata) != 0:
            keys = ['model','error', 'comment']
            #use fitter object for this model to get fit param names
            X = Fitting.makeFitter(fitdata)
            if X != None:
                mdata = X.getFitDict()
            else:
                mdata = {}
            for k in keys:
                if fitdata.has_key(k):
                    mdata[k] = fitdata[k]

        if hasattr(self, '__meta_data__'):
            if self.__meta_data__.has_key(dataset):
                for k in self.__meta_data__[dataset].keys():
                    mdata[k] = self.__meta_data__[dataset][k]
        return mdata

    def getMeta(self, dataset, field):
        """Get a field from the meta data"""
        if self.__meta_data__.has_key(dataset) and self.__meta_data__[dataset].has_key(field):
            return self.__meta_data__[dataset][field]
        else:
            return None

    def setMetaData(self, meta):
        self.__meta_data__ = meta
        return

    def printMeta(self):
        """print current meta data"""
        if not hasattr(self, '__meta_data__'):
            print 'No Meta data'
            return
        M = MetaData(data=self.__meta_data__)
        M.prettyPrint()
        return

    def getField(self, dataset, field):
        """Some fields might be in the ekin dataset dict itself"""
        #we should really use meta_data instead in future!
        if self.data[dataset].has_key(field):
            return self.data[dataset][field]
        else:
            return None

    def fitDataset(self, dataset, model=None, update=True,
                    noiter=300, conv=None, grad=None, silent=False,
                    guess=True, callback=None):
        """Calculate the fit for this current dataset, if a model
           is given we use that instead of the current one.
           update=True means that the dataset fit info will be overwritten"""

        datatofit = self.data[dataset]
        if model == None:
            currfitdata = self.getFitData(dataset)
            if currfitdata.has_key('model'):
                model = currfitdata['model']
            else:
                model = self.defaultmodel
        else:
            currfitdata = None
        fitresult, X = Fitting.doFit(datatofit, fitdata=currfitdata, model=model,
                                     noiter=noiter, conv=conv, grad=grad, silent=silent,
                                     guess=guess, callback=callback)
               
        if fitresult == None:
            print 'Fitter returned None..'
            return

        if update == True:
            self.__datatabs_fits__[dataset] = fitresult

        return fitresult

    def fitDatasets(self, datasets='ALL', update=True, findbest=False, models=None,
                        noiter=300, conv=1e-6, grad=1e-6, silent=False):
        """Fit multiple datasets in this project, if no datasets are given, we
           fit all of them.
           If findbest=True, finds the best fit model using f-testing for the given
           dataset.
           models is optional and provides a list of models to try for best fitting
        """
        datasets = self.checkDatasetParam(datasets)
        print 'fitting %s datasets..' %len(datasets)
        if findbest == True:
            for name in datasets:
                self.findBestModel(name, update=update, models=models)
        else:
            for name in datasets:
                self.fitDataset(name, update=update, model=models[0],
                                noiter=noiter, conv=conv, grad=grad,
                                silent=silent)
        return

    def findBestModel(self, dataset, update=True, models=None, mode=None,
                            checkfunc=None, conv=None, grad=None, alpha=0.01,
                            silent=False):
        """Finds the best fit model using f-testing"""
      
        thisdata = self.data[dataset]
        print 'models', models
        if models == None:
            import copy
            models = copy.deepcopy(self.mode_definition[self.currentmode])

        fitdata, p = Fitting.findBestModel(thisdata, models, checkfunc=checkfunc, conv=conv,
                                            grad=grad, alpha=alpha, silent=silent)

        if update == True:
            if fitdata == None:
                self.setFitData(dataset)
            else:
                self.setFitData(dataset, fitdata)
        return fitdata, p

    def estimateExpUncertainty(self, dataset, xuncert=0.01, yuncert=0.01,
                                runs=10, doplots=False, conv=1e-6, grad=1e-6,
                                guess=True):
        """Chresten's method for estimating the error on each parameter due
           to the effect of the exp uncertainty, which is user provided """

        data = self.data[dataset]
        fitdata = self.__datatabs_fits__[dataset]
        #print 'estimating exp uncertainty for %s' %dataset
        fitstats = Fitting.estimateExpUncertainty(data, fitdata,
                                        xuncert=xuncert, yuncert=yuncert,
                                        runs=runs, doplots=doplots,
                                        conv=conv, grad=grad, guess=guess)

        return fitstats

    def setPylab(self, size=(6,4)):
        """Setup plotting variables if pylab available"""
        try:
            fig = plt.figure()
        except:
            return False
        return True


    def doPlot(self, ek, fitdata, options):
        """Need a bare bones plotter that can plot any given ekin data
           Should then call this one in plotDatasets below and simplify"""
        x,y,a, xerr,yerr = ek.getAll()
        line = ax.scatter(x, y, marker='o')
        #plotting fit should also be seperate?
        f = fitdata
        d=[]
        X = Fitting.makeFitter(f, zip(x,y))
        if X==None: return
        #we now use this model to draw the fit line
        fitclr = lineclr
        if plotoption != 3:
            fitclr='r'
        xx, fity = self.updateFit(X, ax, x, clr=fitclr, plotoption=plotoption)
        return

    def plotDatasets(self, datasets='ALL', data=None, fitdata=None,
                           filename=None, plotoption=1, cols=0,
                           size=(6,4), linecolor=None, figure=None,
                           showfitvars=False, dpi=80, **kwargs):
        """Plot a dataset or list of datasets, if none given all are  plotted.
           plotoptions:
           1 - each dataset in one plot, different figures/images
           2 - each dataset in different axes but same figure - one image
           3 - all datasets in one plot, single figure
           Note: errorbars are drawn +/- value supplied"""
        status = self.setPylab(size)
        prms = Params(**kwargs)
        self.current = datasets        
        shapes = Options.shapes

        plt.rc('font', family=prms.font)
        plt.rc('font', size=prms.fontsize)        
        plt.rc('text', usetex=prms.usetex)

        if plotoption == 2 and cols>=3:
            size=(3,4)
        if status == False:
            print 'sorry, install pylab..'
            return
        if figure != None:
            fig = figure
        else:
            fig = plt.figure(figsize=size, dpi=80)
        fig.clf()
        datasets = self.checkDatasetParam(datasets)
        if plotoption == 2:
            if cols==0 or cols=='':
                cols=math.ceil(math.sqrt(len(datasets)))               
            plt.rc('font', size=9.0)
            if prms.title != None:
                fig.suptitle(prms.title, fontsize=prms.fontsize)
            if len(datasets) == 2:
                   cols=1; dim=2
            else:       
                dim = math.ceil(len(datasets)/float(cols))
            n=1            
        else:
            ax = fig.add_subplot(111)
            ax.cla()
            self.ax = ax

        legendlines = []
        a = numpy.arange(-4,4,0.5).tolist()
        #get ekin data into plot values
        xdata={};ydata={};adata={};xerrors={};yerrors={}

        for name in datasets: 
            self.fitline = None
            ek = self.data[name]            
            xd,yd = xdata[name],ydata[name] = ek.getxy()
            adata[name] = ek.active
            xerrors[name], yerrors[name] = ek.getErrors()
            #if no error stored for points, try using global error if present
            if xerrors[name] == None and xerror != None:
                xerrors[name] = []
                for i in range(len(xd)): xerrors[name].append(xerror)
            if yerrors[name] == None and yerror != None:
                yerrors[name] = []
                for i in range(len(yd)): yerrors[name].append(yerror)

            if prms.xlabel != '':
                xlabel = prms.xlabel 
            else:
                xlabel = ek.labels[0]            
            if prms.ylabel != '':
                ylabel = prms.ylabel 
            else:
                ylabel = ek.labels[1]  
 
        e=[]
        for i in xdata: e.append(i[0])
        i=0
        cc=0
        for name in datasets:
            if cc >= len(prms.colors):
                cc=0
            if prms.varyshapes == True:
                prms.marker = shapes[cc]
            if plotoption == 2:
                ax = fig.add_subplot(int(dim),int(cols),n)
                ax.set_title(name, fontsize=prms.fontsize, fontweight='bold')
                if cols<4:
                    ax.set_xlabel(xlabel,family=prms.font)
                    ax.set_ylabel(ylabel,family=prms.font)
                if prms.grid == True:
                    ax.grid(True)
                n=n+1
                self.ax = ax

            x = xdata[name]
            y = ydata[name]
            act = adata[name]
            xerr = xerrors[name]
            yerr = yerrors[name]
            self.currxdata = x
            if len(x) == 0:                
                continue
            if prms.normalise == True:
                self.mindata=min(y); self.maxdata=max(y)
                y = normList(y, inputmin=self.mindata, inputmax=self.maxdata)
            if prms.graphtype == 'XY':
                if prms.logx == True and prms.logy == False:
                    line, = ax.semilogx(x, y, prms.marker, alpha=prms.alpha)
                elif prms.logy == True and prms.logx == False:
                    line, = ax.semilogy(x, y, prms.marker, alpha=prms.alpha)
                elif prms.logx == True and prms.logy == True:
                    line, = ax.loglog(x, y, prms.marker, alpha=prms.alpha)
                elif prms.marker in ['-','--',':','.',':.',',','|','*']:
                    line, = ax.plot(x, y, prms.marker, linewidth=prms.linewidth, mew=prms.linewidth,
                                        alpha=prms.alpha, ms=prms.markersize)
                else:
                    ptcolors=[]
                    if prms.grayscale==True:
                        actclr = '0.5'
                    else:
                        actclr = prms.colors[cc]
                    for i in range(len(x)):              
                        if act[i]==0:
                            ptcolors.append('0.8')
                        else:
                            ptcolors.append(actclr)
                          
                    line = ax.scatter(x, y, marker=prms.marker, c=ptcolors,
                                      s=prms.markersize, lw=prms.linewidth, alpha=prms.alpha)
                    cc+=1

                lineclr = self._getlineColor(line)
                if prms.showerrorbars == True and (yerr != None or xerr != None):
                    #print xerr,yerr
                    errline = ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=None,
                                            elinewidth=1.0, ecolor=lineclr, alpha=prms.alpha)
            elif prms.graphtype == 'bar':
                bars = ax.bar(x, y, linewidth=prms.linewidth, alpha=prms.alpha)
                line = bars[0]
                lineclr = line.get_facecolor()

            if len(xd)>0:
                if abs(min(yd))>0.01 and max(yd)<1e4:               
                    major_formatter = plt.FormatStrFormatter('%2.2f')
                    ax.yaxis.set_major_formatter(major_formatter)

            legendlines.append(line)
            if self.__datatabs_fits__.has_key(name):
                f = self.__datatabs_fits__[name]
                if f.has_key('model'):                    
                    #create the required fitter
                    d=[]
                    X = Fitting.makeFitter(f, zip(x,y))
                    if X==None:
                        continue
                    #we now use this model to draw the fit line
                    fitclr = lineclr
                    if plotoption != 3:
                        fitclr='r'
                    if prms.grayscale==True:
                        fitclr='0.4'
                    xx, fity = self.updateFit(X, ax, xdata[name], clr=fitclr,
                                              normalise=prms.normalise,
                                              plotoption=plotoption)

            if prms.xlim!=0:
                l = min(x)-prms.xlim; u=max(y)+prms.xlim
                ax.set_xlim((l,u))
                print l,u
            if prms.ylim!=0: 
                l = min(y)-prms.ylim; u=max(y)+prms.ylim
                ax.set_ylim((l,u))
                               
        if plotoption != 2:
            if prms.title == None or prms.title == '':
                prms.title=name
            ax.set_title(prms.title,family=prms.font)
            ax.set_xlabel(xlabel,family=prms.font)
            ax.set_ylabel(ylabel,family=prms.font)

            if showfitvars == True:
                self.showfitResults(X, ax)
            if prms.grid == True:
                ax.grid(True)
            if prms.legend == True:
                l=ax.legend(legendlines, datasets, numpoints=1,
                            loc=prms.legendloc, prop=FontProperties(size="smaller"))
                l.get_frame().set_alpha(0.8)
        else:
            fig.subplots_adjust(hspace=0.4,wspace=0.4)

        if filename != None:
            fig.savefig(filename, dpi=dpi)
            plt.close(fig)

        return ax
        
    def updateFit(self, X, ax=None, xdata=None, clr=None,
                        normalise=None, plotoption=1):
        """Get the fit data line for given fitter object and redraw"""
        if xdata == None:
            xdata = self.currxdata
        if ax == None:
            ax = self.ax
        if clr == None:
            clr='r'
        e=[]
        for i in xdata: e.append(i)
        fity=[]
        inc = (max(e)-min(e))/50
        xx = numpy.arange(min(e)-inc, max(e)+inc, inc).tolist()
        for i in xx:
            fity.append(X.evaluate(i))
        if normalise == True:
            fity = normList(fity, inputmin=self.mindata, inputmax=self.maxdata)
        if clr==None:
            clr=line.get_color()
        if self.fitline == None or plotoption != 1:
            self.fitline, = ax.plot(xx, fity, clr, linewidth=2.5, alpha=0.8)
        else:
            self.fitline.set_data(xx, fity)
        return xx, fity

    def _getlineColor(self, line):
        """Private method - get a lines' color"""
        clr = '#cccccc'
        try:
            clr = line.get_color()
        except:
            i=0
            while clr == '#cccccc':
                clr = tuple(line.get_facecolor()[i][:3])
                clr = matplotlib.colors.rgb2hex(clr)
                i+=1        
        return clr
    
    def showfitResults(self, X, ax):
        """Show fit vars in a plot"""

        ax.text(.1, .8, X.getResult(), transform=ax.transAxes,
                bbox={'facecolor':'yellow', 'alpha':0.6, 'pad':10})
        plt.draw()
        return

    def updatePlot(self):
        """Update current plot points"""
        curr = self.current
        if type(curr) is not types.StringType:
            return
        #data =
        l = self.ax.get_lines()
        #print l
        return

    def checkDatasetParam(self, datasets):
        """Check dataset parameter to see if it's a list or not, allows
           methods to be called with dataset as a string or list"""
        if datasets=='ALL':
            return self.datasets
        else:
            import types
            if type(datasets) is types.StringType:
                s=datasets
                datasets=[]
                datasets.append(s)
                return datasets
            elif type(datasets) is types.ListType:
                return datasets
            else:
                print 'unknown format for dataset names, provide a list'
                return None

    def __repr__(self):
        """Overload string repr of class"""
        l = len(self.datasets)
        return "Ekin project containing %s datasets" %l

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, item):
        self.__dict__[key] = item
        
    def getFields(self):
        return self.__dict__.keys()
    
    def addProject(self, E, overwrite=False):
        """Combine another project with the current one, duplicates are renamed"""
        import copy
        for d in E.datasets:            
            edata = copy.deepcopy(E.getDataset(d))            
            fitdata = E.getFitData(d)
            if d in self.datasets:
                name = d + '_1'
            else:
                name = d            
            self.insertDataset(edata, name, fit=fitdata)
            for field in self.ekinfields:
                if not E.__dict__.has_key(field):
                    continue
                #if E.__dict__.has_key(field):
                #    print field, type(E.__dict__[field]), E.__dict__[field][name]                
                if E.__dict__[field].has_key(name):   
                    self.__dict__[field][name] = copy.deepcopy(E.__dict__[field][name])
        return

#
# These methods are dataset based and could be handled by the ekindataset class
#
    def setAxesLabels(self, x, y):
        """Set axes labels for all datasets"""
        for name in self.datasets:
            self.data[name][0]['label'] = x
            self.data[name][1]['label'] = y
        return

    def getAxesLabels(self, dataset):
        """Get axes labels"""
        try:
            labels = (self.data[dataset][0]['label'], self.data[dataset][1]['label'])
            return labels
        except:
            return None

class Params(object):
    """Class to store options passed to plotter"""
    def __init__(self, **kwargs):        
        self.colors = Options.colors
        self.marker = 'o'
        self.linewidth = 1
        self.legend = False
        self.legendloc = 'best'
        self.font = 'sans-serif'
        self.fontsize = 12
        self.graphtype = 'XY'
        self.grid = False
        self.alpha = 0.8
        self.usetex = False
        self.markersize = 25
        self.normalise = False
        self.showerrorbars = False
        self.logx = False
        self.logy = False
        self.title = None
        self.xlabel = ''
        self.ylabel = ''
        self.varyshapes = False
        self.grayscale = False        
        for k in kwargs:
            self.__dict__[k] = kwargs[k]
        return   

