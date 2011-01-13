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
from Tkinter import *
import Pmw
import os, sys, math, random, glob, numpy, string
import ConfigParser, csv, xlrd
from PEATDB.Ekin.Base import EkinProject, EkinDataset
from PEATDB.Ekin.Web import EkinWeb
#from PEATDB.Ekin.Convert import EkinConvert
from PEATDB.Ekin.Fitting import Fitting
from PEATDB.Optparse import OptparseDialog

class KineticsAnalyser(Plugin):
    """A class for processing specific Novozymes Kinetics exp data"""
    """Author: Damien Farrell"""

    capabilities = ['gui','uses_sidepane']
    requires = ['numpy','matplotlib']
    menuentry = 'Kinetics Analysis'
    gui_methods = []

    times=range(0,910,10)
    abserr = 0.005      #absorbance perc error
    pherr = 0.1        #ph absolute error
    sconcerr = 5     #substrate perc conc error
    inhibconcerr = 1 #inhib conc perc error
    sets = ['B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']
    temps = {4:'room',5:'60',6:'65',7:'70',8:'75'}
    bold = "\033[1m"; reset = "\033[0;0m"

    def __init__(self, parent=None, db=None,
                        astdef=None, mmdef=None):
        if parent!=None:
            self.parent = parent
        self.DB = db
        if mmdef != None:
            self.parseMMdefs(mmdef)
        self.parser = self.createOptions()
        return

    def main(self, parent):
        if parent==None:
            return
        self.parent = parent
        self.DB = parent.DB
        self.parentframe = None
        self._doFrame()
        return

    def _doFrame(self):
        self.ID='Kinetics Data Analysis Plugin'
        if 'uses_sidepane' in self.capabilities:
            self.mainwin = self.parent.createChildFrame()
        else:
            self.mainwin=Toplevel()
            self.mainwin.title(self.ID)
            self.mainwin.geometry('800x600+200+100')

        methods = self._getmethods()
        methods = [m for m in methods if m[0] in self.gui_methods]
        l=Label(self.mainwin, text='Kinetics/Temp stab data analysis')
        l.pack(side=TOP,fill=BOTH)
        b=Button(self.mainwin, text='Execute',
               command= self.callfromGUI, bg='#FFFF66')
        b.pack(side=TOP,fill=BOTH)
        b=Button(self.mainwin, text='Close', command= self.close)
        b.pack(side=TOP,fill=BOTH)
        optsdlg = OptparseDialog()
        win, self.optvars = optsdlg.createWidgets(self.mainwin, self.parser)
        win.pack(side=LEFT,fill=BOTH)
        from textFrame import HyperlinkManager
        self.log = Pmw.ScrolledText(self.mainwin,
                labelpos = 'n',
                label_text='Logs',
                usehullsize = 1,
                hull_width = 800,
                hull_height = 500,
                text_wrap='word')
        self.log.pack(side=LEFT,fill=BOTH,padx=4,pady=4)
        hyperlink = HyperlinkManager(self.log)
        self.log.insert(INSERT, "see help here", hyperlink.add_link('http://enzyme.ucd.ie'))
        return

    def callfromGUI(self):
        """Run from GUI, passing opt vars"""
        opts, remainder = self.parser.parse_args()
        for o in opts.__dict__.keys():
            if o in self.optvars:
                opts.__dict__[o] = self.optvars[o].get()

        #redirect stdout to log win
        self.log.delete(1.0,END)
        sys.stdout = self
        sys.stderr = self
        self.run(opts)
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        return

    def write(self, txt):
        """Handle stdout"""
        self.log.insert(END, txt)
        self.log.update_idletasks()
        return

    def close(self):
        self.mainwin.destroy()
        return


# ---------------

    def doHeader(self):
        print '<html>'
        print '<head>'
        print '<script type="text/javascript" src="sortable.js"></script>'
        print '<link href="../../table.css" rel="stylesheet" type="text/css">'
        print '</head>'
        print '<body>'
        return

    def parseDefs(self, filename):
        """Parse a basic def file"""
        self.defsfilepath = os.path.dirname(filename)
        f = open(filename,'r')
        c = ConfigParser.ConfigParser()
        c.read(filename)
        fields = ['path', 'setname','phlist','subconcs','inhibconcs','variantnames','required',
                    'runs', 'astdate', 'astdilfile', 'rawdilfile', 'otherdilutions']

        for field in fields:
            try:
                self.__dict__[field]=eval(c.get('default',field))
            except:
                self.__dict__[field]=str(c.get('default',field))

        if self.required == '':
            self.required = self.variantnames
        print 'parsed defs file, got these params:'
        print '----------------------------------'
        for f in fields:
            print self.bold + f + self.reset, self.__dict__[f]
        print '----------------------------------'
        return


    def doplots(self, filename=None):
        """make plots of fits"""
        #print 'opening %s' %os.path.join(path, filename+'.ekinprj')
        ew = EkinWeb()
        if not os.path.exists(os.path.join(path, 'plots')):
            os.mkdir(os.path.join(path, 'plots'))
        ew.showEkinPlots(project=os.path.join(path, filename+'.ekinprj'),
                                outfile=os.path.join(path, filename+'plots.html'),
                                imgpath=os.path.join(path, 'plots/'),
                                path='plots/',
                                title=filename+' fits',
                                columns=3)


    def getDilutions(self, filename, xls=False):
        """get dilutions from txt or xls file"""
        dilutions = {}
        dilutions2 = {}
        line=0
        self.template_enzconcs={}
        if xls == False:
            reader = csv.reader(open(filename), dialect='excel-tab')
            for row in reader:
                if len(row) == 0 or len(row[1]) ==0 : continue
                name = row[0]
                if line<=11:
                    dilutions[name] = float(row[1])
                    if len(row)>2:
                        self.template_enzconcs[name]=float(row[2])
                else:
                    dilutions2[name] = float(row[1])
                line+=1
        else:
            book = xlrd.open_workbook(filename)
            sh1 = book.sheet_by_index(0)
            for r in range(12):
                name=sh1.col_values(0)[r]
                dilutions[name] =  sh1.col_values(1)[r]
        print 'got dilutions from %s' %filename
        print 'dilutions1:', dilutions
        print 'dilutions2:', dilutions2
        return dilutions, dilutions2


    def loadTempData(self):
        """Import temp data from 2 xls files for each set
         set1.xls contains duplicate data for variants 1-6
         i.e. columns 1-6 and 7-12 are repeats for 6 variants"""

        print 'importing temp data..'

        filename1 = os.path.join('temp', 'set1.xls') #variants 1-6
        filename2 = os.path.join('temp', 'set2.xls') #variants 7-12
        book1 = xlrd.open_workbook(os.path.join(self.path, filename1))
        try:
            book2 = xlrd.open_workbook(os.path.join(self.path, filename2))
        except:
            book2 = None


        rawdata={}
        times=range(0,190,10)
        mins=[]
        for t in times: mins.append(round(t/60.0,3))
        rowoffset=9
        i=0

        for i in range(0,12):
            name = self.variantnames[i]
            rawdata[name] = {}
            for dupname in ('1','2'):
                rawdata[name][dupname] = {}

        coloffset=0
        for book in (book1, book2):
            if book == None: continue
            for sheet in range(4,9):
                i=coloffset
                sh = book.sheet_by_index(sheet)
                temp = self.temps[sheet]
                for dup in (0,6): #7-12 are duplicates
                    i=i-dup
                    if dup == 0: dupname = '1'
                    else: dupname = '2'
                    #get raw data from cols 1-6 in sheet
                    for c in range(2+dup,8+dup):
                        row=10
                        name = self.variantnames[i]
                        #print book, sh, dup, c, name
                        rawdata[name][dupname][temp] = {}
                        for s in range(0,8):
                            ph = self.phlist[s]
                            r=row
                            yvals=[]; yerrs=[]; xerrs=[]
                            for t in times:
                                #print s, t, r, c, sh.col_values(c)[r], temp, ph, name
                                val=sh.col_values(c)[r]
                                yvals.append(val)
                                try:
                                    yerrs.append(val*self.abserr)
                                except:
                                     yerrs.append(0.0)
                                xerrs.append(0.1)
                                r = r + rowoffset
                            row+=1
                            rawdata[name][dupname][temp][ph] = EkinDataset(xy=[mins, yvals],
                                                errorlists=(xerrs, yerrs), labels=('t','A'))
                        i+=1
            coloffset+=6

        '''for n in rawdata.keys():
            for d in rawdata[n].keys():
                print n, d, rawdata[n][d].keys()'''

        return rawdata

    def processTempData(self, raw):
        """Process raw temperature data"""
        print 'processing temp data..'
        import copy
        actlists={}

        for name in self.required:
            print name
            E = EkinProject(mode='General')
            for d in raw[name].keys():
                for t in raw[name][d].keys():
                    for p in raw[name][d][t].keys():
                        E.insertDataset(raw[name][d][t][p], d+'_'+t+'_'+str(p))
            E.fitDatasets('ALL', models=['Linear'], noiter=100, conv=1e-6, grad=1e-8, silent=True)
            for d in E.datasets:
                status = self.checkLinearFit(E, d, cutoff = 0.005, fromstart=False, percs=[0.7, 0.5])
                actlists[d] = status
            print 'estimating exp uncertainties..'
            for d in E.datasets:
                ferrs = E.estimateExpUncertainty(d, runs=10)
                E.addMeta(d, 'exp_errors', ferrs)
            E.saveProject(os.path.join(self.path, name.replace(' ','')+'_temp'))

        tempslst = ['room','60','65','70','75']
        for name in self.required:
            #we create a plot for each variant of vel(slope) vs. temp @ each pH
            E = EkinProject()
            E.openProject(os.path.join(self.path, name.replace(' ','')+'_temp'))
            datasets = copy.deepcopy(E.datasets)
            if len(datasets) == 0:
                continue

            for dup in ('1','2'):
                phvals = []; tms = []; tmerrs = []
                for ph in self.phlist:
                    active = []
                    temps = [];activities=[];xerrors=[];yerrors=[]
                    ref = dup+'_'+'room'+'_'+ph
                    refact = E.getMetaData(ref)['a']
                    refacterr =  E.getMeta(ref, 'exp_errors')['a'][1]
                    for temp in tempslst:
                        d = dup+'_'+str(temp)+'_'+ph

                        act = E.getMetaData(d)['a']
                        if act < 0:
                            print 'negative slope'
                            act=0
                        acterr = E.getMeta(d, 'exp_errors')['a'][1]
                        percact = act / refact
                        yerror = percact * math.sqrt( pow(acterr/act,2) + pow(refacterr/refact,2))
                      
                        active.append(actlists[d])
                        if temp == 'room':
                            temp='25'

                        temps.append(temp)
                        activities.append(percact)
                        yerrors.append(yerror)
                    if activities[1]>1.0:
                        print ph, refact, activities[1]
                        active[0]=0
                    #print name, ph, active
                    dsetname = dup+'_'+ph
                    ek = EkinDataset(xy=[temps, activities], yerrs=yerrors, active=active,
                                                xlabel='temp',ylabel='res act.')
                    E.insertDataset(ek, dsetname)
                    E.fitDataset(dsetname, model='Residual Activity', conv=1e-6, silent=True)
                    try:
                        ferrs = E.estimateExpUncertainty(dsetname, runs=10)
                        E.addMeta(dsetname, 'exp_errors', ferrs)
                    except:
                        ferrs = E.estimateExpUncertainty(dsetname, runs=2)
                        E.addMeta(dsetname, 'exp_errors', ferrs)
                    phvals.append(float(ph.replace(',','.')))
                    tms.append(E.getMetaData(dsetname)['t50'])
                    tmerrs.append(E.getMeta(dsetname, 'exp_errors')['t50'][1])

                #finally add temp vs ph plot and fit
                ek = EkinDataset(xy=[phvals, tms], xerr=None, yerr=tmerrs,
                                                xlabel='ph',ylabel='tm')
                E.insertDataset(ek, dup+'_'+'pHvsTm')
                E.fitDataset(dup+'_'+'pHvsTm', model='1 pKa 2 Chemical shifts', conv=1e-6, silent=True)
            E.saveProject()
        return


    def loadRawData(self, path, variantnames, phvalues, date, setname,
                         run=None, subconcs=None, dilutions=None):
        """Load raw MM data for one set of replicates and return dict of variants.
           raw data is loaded from text files containing sets per pH and
           the other args provide the exp conditions
           parseMMdefs uses MM.def to get exp details and provide them here
           """
        print 'processing raw data..'
        rawdata = {}

        for name in variantnames:
            rawdata[name] = {}
            for ph in phvalues:
                rawdata[name][ph]={}

        for ph in phvalues:
            print 'processing raw data for ph', ph

            filename = 'set%s_%s_%s_pH%s.txt' % (setname, date, run, ph)
            filename = os.path.join('MM', filename)
            try:
                reader = csv.reader(open(os.path.join(path, filename)),
                                    dialect='excel-tab')
            except:
                print 'failed to open %s !' %filename
                print
                continue

            rows=[]
            for row in reader:
                rows.append(row)
            mins=[]
            for t in self.times: mins.append(round(t/60.0,3))
            i=0
            row=3
            offset=9
            #get raw absorbance data from text file instead of excel
            for c in range(2,14):
                row=3
                name = variantnames[i]
                for s in range(0,8):
                    r=row
                    yvals=[];yerrs=[]
                    for t in self.times:
                        #print s, t, r, c, len(rows[r]), rows[r][c]
                        try:
                            val=float(rows[r][c].replace(",", "."))
                            yvals.append(val)
                            yerrs.append(val*self.abserr/100)
                        except:
                            yvals.append('')
                            yerrs.append(0)
                            #print 'Nan, skipping'
                        r = r + offset
                    row+=1
                    S=subconcs[s]
                    rawdata[name][ph][S] = EkinDataset(xy=[mins, yvals],
                                               yerrs=yerrs, xlabel='t',ylabel='A')
                i+=1

        return rawdata


    def processKineticsRaw(self, variantdata, phvalues):
        """Get raw data for a variant into ekin projects, fit MM and save"""

        E = EkinProject(mode='General')
        for ph in variantdata.keys():
            for s in variantdata[ph].keys():
                E.insertDataset(variantdata[ph][s], ph+'_'+str(s))

         #inital fit of raw data to get slopes
        print 'fitting raw velocity data..'
        E.fitDatasets(models=['Linear'], noiter=100, conv=1e-6, grad=1e-8, silent=True)

        #check linear fits applying criteria
        print 'checking linear fits..'
        actlists={}
        for ph in phvalues:
            actlists[ph]=[]
            for d in E.datasets:
                x = d.split('_')
                p, sconc = x[0], round(float(x[1]),3)
                if p == ph:
                    status = self.checkLinearFit(E, d, cutoff = 0.04)
                    actlists[ph].append(status)
        print actlists
        #estimate exp uncertainty and save linear fit errors in the meta_data
        print 'estimating exp uncertainties for raw data..'
        for d in E.datasets:
            ferrs = E.estimateExpUncertainty(d, runs=10)
            E.addMeta(d, 'exp_errors', ferrs)

        #Extract mm data and fit to michealis-menten
        print 'creating MM data..'
        import copy
        datasets = copy.deepcopy(E.datasets)
        mmdatasets = []
        for ph in phvalues:
            slopes=[];sconcs=[];xerrors=[];yerrors=[]
            active = actlists[ph]
            for d in datasets:
                x = d.split('_')
                p, sconc = x[0], round(float(x[1]),3)
                #print ph, name, n, sconc, slope, yerror
                if p == ph:
                    #get slope and error and add the point
                    try:
                        slope = E.getMetaData(d)['a']
                        yerror = E.getMeta(d, 'exp_errors')['a'][1]
                    except:
                        print 'no slope, skipping'
                        continue
                    slopes.append(slope)
                    sconcs.append(sconc)
                    xerrors.append(sconc*self.sconcerr/100)
                    yerrors.append(yerror)
            #print ph, name, sconcs, slopes, yerrors
            edata = EkinDataset(xy=[sconcs, slopes], xerrs=xerrors, yerrs=yerrors,
                                        active=active, xlabel='[S]',ylabels='V')

            E.insertDataset(edata, ph)
            mmdatasets.append(ph)

        E.fitDatasets(datasets=mmdatasets, models=['Michaelis-Menten'], noiter=400,
                            conv=1e-7, grad=1e-8, silent=True)

        print 'estimating exp uncertainties for MM data..'
        for d in mmdatasets:
            ferrs = E.estimateExpUncertainty(d, runs=10)
            E.addMeta(d, 'exp_errors', ferrs)

        return E


    def analyseKinetics(self, name, enzymeconc, enzymeconcerr,
                              dilution, altdilution, doftest=True, run=None):
        """Analyse MM kinetics fits from ekin project to extract kcats and kms etc"""

        print 'analysing MM data..'

        E = EkinProject()
        E.openProject(os.path.join(self.path, name))

        #We iterate over the phs to get kcats vs pH etc.
        kcatdatasets = []
        phs=[]; kcats=[]; kms=[]; kcatkms=[]; kslopes=[]
        kcaterrors=[]; kmerrors=[]; kcatkmerrors=[]; pherrors=[]; kslopeerrors=[]
        activelist=[]
        enzpercerr = enzymeconcerr/enzymeconc
        #introduce an error factor for the enzyme conc here to account
        #for variance in kcats between replicates..
        factor = 0.1
        lbls=['ph', 'vmax','[E]', '[E]err', 'dilution', 'kcat', 'kcat err', 'km', 'km err',
                 'kcat/km', 'kcat/km err', 'slope@0.025', 'slope@0.05']

        #send output to log files as html tables
        saveout = sys.stdout
        fsock = open(os.path.join(self.path,'log_'+str(run)+'.html'), 'a')
        sys.stdout = fsock
        self.doHeader()

        print '<h1>%s</h1>' %name
        print '<table>'
        for l in lbls:
            print '<th>%s</th>' %l
        print '<tr>'
        for d in self.phlist:
            ph=float(d.replace(',','.'))
            phs.append(ph)
            if d == '6' or d=='6,5':
                dil = dilution
            else:
                dil = altdilution

            vmax = E.getMetaData(d)['Vmax']
            km = E.getMetaData(d)['Km']
            kms.append(km)

            kcat = vmax / ( enzymeconc / dil )
            kcatkm = (kcat/km) / 1000

            vmaxerror = E.getMeta(d, 'exp_errors')['Vmax'][1]
            kmerror = E.getMeta(d, 'exp_errors')['Km'][1]
            #kmpercerror = kmerror / (km / dil)

            #check MM fit
            status = self.checkFit(E, d, 0.002)
            activelist.append(status)

            #get slope at chosen [S] to use for alternate kcat/km value
            #we check that it's less than km at this [S]
            ekd = E.getDataset(d)
            kslope1 = ekd.getYVal(4) #0.025
            kslope2 = ekd.getYVal(5) #0.05

            kslopes.append(kslope1)
            kslopeerrors.append(ekd.getYError(4))

            kcaterror = kcat * math.sqrt( pow(vmaxerror/vmax,2) + pow(enzpercerr+factor,2) )
            kcatkmerror = kcatkm * math.sqrt( pow(kcaterror/kcat,2) + pow(kmerror/km,2) )
            kcats.append(kcat)
            kcaterrors.append(kcaterror)
            kmerrors.append(kmerror)
            pherrors.append(self.pherr)
            kcatkms.append(kcatkm)
            kcatkmerrors.append(kcatkmerror)

            prms=[ph, vmax, enzymeconc, enzymeconcerr, dil, kcat, kcaterror, km, kmerror, kcatkm, kcatkmerror]
            for p in range(len(prms)):
                print '<td>%.3f</td>' %float(prms[p])
            #print 'using slope %s at [S]=%s' %(kslope, ekd.getXVal(5))
            print '<td>%.4f</td>' %float(kslope1)
            print '<td>%.4f</td>' %float(kslope2)
            print '<tr>'

        print '</table>'

        sys.stdout.flush()
        sys.stdout = saveout
        fsock.close()

        #print activelist
        kcatdata = EkinDataset(xy=[phs, kcats], xerrs=pherrors, yerrs=kcaterrors,
                                        active=activelist,xlabels='pH',ylabels='Kcat')

        kmdata = EkinConvert.xy2ekin([phs, kms], errorlists=(pherrors, kmerrors), labels=('pH','Km'))
        kcatkmdata = EkinConvert.xy2ekin([phs, kcatkms], errorlists=(pherrors, kcatkmerrors), labels=('pH','Kcat/Km'))
        kslopesdata = EkinConvert.xy2ekin([phs, kslopes], errorlists=(pherrors, kslopeerrors), labels=('pH','Kcat/Km'))
        E.insertDataset(kcatdata, 'Kcat', replace=True)
        E.insertDataset(kmdata, 'Km',replace=True)
        E.insertDataset(kcatkmdata, 'Kcat/Km',replace=True)
        E.insertDataset(kslopesdata, 'Slope@0.025',replace=True)
        E.addMeta('Kcat', 'enzymeconc', enzymeconc)

        kd = ['Kcat', 'Km', 'Kcat/Km', 'Slope@0.025']

        if doftest == True:
            #use ftest to see if 1 pka is better than fitting a line
            for d in kd:
                E.findBestModel(d, update=True, models=['Linear', '1 pKa 2 Chemical shifts'],
                                         conv=1e-12, grad=1e-12, strictchecking=False)
        else:
            for d in kd:
                E.fitDatasets(datasets=['Kcat', 'Km', 'Kcat/Km', 'Slope@0.025'], models=['1 pKa 2 Chemical shifts'],
                                    noiter=400, conv=1e-12, grad=1e-12, silent=True)
        print 'got pH vs. Kcat data'
        print 'done.', E
        print
        return E

    def combineKcatProjects(self):
        """put kcats, kms and kcats/kms for all variants into single projects"""
        print 'creating kcat/km projects..'
        Ekcats = EkinProject()
        Ekms = EkinProject()
        Ekcatskms = EkinProject()
        Ekslopes = EkinProject()
        required = self.required
        runs = self.runs
        kd = {'Kcat':Ekcats, 'Km':Ekms, 'Kcat/Km':Ekcatskms, 'Slope@0.025':Ekslopes}
        for name in required:
            for r in runs:
                E = EkinProject()
                E.openProject(os.path.join(self.path, name.replace(' ','')+'_'+str(r)))
                for k in kd:
                    edata = E.getDataset(k)
                    fdata = E.getFitData(k)
                    kd[k].insertDataset(edata, name+'_'+str(r))
                    kd[k].setFitData(name+'_'+str(r), fdata)

                #estimating exp error for fit in kcats
                ferrs= Ekcats.estimateExpUncertainty(name+'_'+str(r), runs=5)
                Ekcats.addMeta(name+'_'+str(r), 'exp_errors', ferrs)

        Ekcats.saveProject(os.path.join(self.path, 'kcats'))
        Ekms.saveProject(os.path.join(self.path, 'kms'))
        Ekcatskms.saveProject(os.path.join(self.path, 'kcats_kms'))
        Ekslopes.saveProject(os.path.join(self.path, 'kslopes'))
        self.Ekcats = Ekcats
        return Ekcats, Ekms, Ekcatskms, Ekslopes


    def loadASTData(self, path, variantnames, date, setname, inhibconcs=None):

        """requires us to get the dilution factor for each conc"""

        filename = 'set%s_%s_ast.xls' % (setname, date)
        filename = os.path.join('AST', filename)
        print filename

        book = xlrd.open_workbook(os.path.join(path, filename))
        sh1 = book.sheet_by_index(0)
        rawdata={}
        times=range(0,310,10)
        mins=[]
        for t in times: mins.append(round(t/60.0,3))
        i=0
        row=3
        offset=9
        #get raw enzyme abs data from sheet and put into ekin
        for c in range(2,14):
            row=3
            name = variantnames[i]
            rawdata[name] = {}
            for s in range(0,8):
                r=row
                yvals=[]; yerrs=[]
                for t in times:
                    #print s, t, r, c, sh1.col_values(c)[r]
                    val=sh1.col_values(c)[r]
                    yvals.append(val)
                    if val == '': print 'val is empty'
                    yerrs.append(val*self.abserr)
                    r = r + offset
                row+=1
                I=inhibconcs[s]
                rawdata[name][I] = EkinConvert.xy2ekin([mins, yvals],
                                    errorlists=(None, yerrs), labels=('t','A'))
            i+=1

        return rawdata


    def fitAST(self, variantdata=None, dilution=None, Eraw=None):
        """Fit velocities vs inhibitor conc. to get enzyme concentrations,
           we get the error on the enz conc from the exp uncert of fit"""

        E = EkinProject(mode='General')
        for i in variantdata.keys():
            E.insertDataset(variantdata[i], str(i))

        #inital fit of raw data to get slopes
        E.fitDatasets(models=['Linear'], noiter=200, conv=1e-8, grad=1e-8, silent=True)


        #check fits applying criteria and deactivate points accordingly
        #criteria: V>0, leave out 3.2 val, check err on slope
        active=[]
        for d in E.datasets:
            status = self.checkLinearFit(E, d, cutoff = 0.08, fromstart=False)
            ferrs = E.estimateExpUncertainty(d, runs=10)
            E.addMeta(d, 'exp_errors', ferrs)
            #print d, ferrs
            slope = E.getMetaData(d)['a']
            ek=EkinDataset(E.getDataset(d))
            Yvals = ek.getY()
            yrange = max(Yvals)-min(Yvals)
            yerror = E.getMeta(d, 'exp_errors')['a'][1] / yrange
            if d == '3.2' or slope < 0.0002 or yerror > 0.02:
                active.append(0)
            else:
                active.append(1)
            #print d, slope, yerror, active
        print active
        print 'creating Inhib conc data..'
        import copy
        datasets = copy.deepcopy(E.datasets)

        #extract slopes and create Inhib data
        slopes=[];iconcs=[]
        xerrors=[];yerrors=[]
        for d in datasets:
            x = d.split('_')
            iconc = round(float(x[0]),3)
            #get slope and error and add the point
            slope = E.getMetaData(d)['a']
            yerror = E.getMeta(d, 'exp_errors')['a'][1]
            slopes.append(slope)
            iconcs.append(iconc)
            xerrors.append(iconc*self.sconcerr/100)
            yerrors.append(yerror)

        #print ph, name, sconcs, slopes, yerrors
        edata = EkinConvert.xy2ekin([iconcs, slopes], errorlists=(xerrors, yerrors),
                                        activelist=active, labels=('[I]','V'))

        E.insertDataset(edata, 'enzyme_conc')
        enzymeconc = self.getEnzymeConcfromAST(E, dilution)

        return E, enzymeconc

    def fitAST2(self, variantdata=None, dilution=None, Eraw=None):
        """Fit AST using Lars method, we fit each raw abs data to MM by first
        transforming Abs vals to [S] and getting a Vmax for each [I]"""
        import copy

        E = EkinProject(mode='General')
        for i in variantdata.keys():
            E.insertDataset(variantdata[i], str(i))

        #transform dataset in terms of [S] vs V vals using Ainit, Amax
        #using eq: [S] = [S]init * (1 - (A405 - A405,init) / (A405,max - A405,init))
        Sinit = 1
        #get Amax?
        Es = EkinProject()

        #get Amax from highest well
        Amax = 0
        for d in E.datasets:
            x,y,a = EkinConvert.ekin2xy(E.getDataset(d))
            if max(y) > Amax:
                Amax = max(y)
        print 'Amax', Amax
        for d in E.datasets:
            x,y,a = EkinConvert.ekin2xy(E.getDataset(d))
            S=[]; V=[]
            Ainit = 0.05
            for a in range(len(y)-1):
                dT=x[a+1]-x[a]
                A = (y[a+1] + y[a])/2
                dA = (y[a+1] - y[a])    #velocity is rate of change Abs

                V.append(dA/dT)
                S.append(abs(Sinit * (1- ((A-Ainit) / (Amax-Ainit)))))

            sdata = EkinConvert.xy2ekin([S,V], labels=('[S]','v'))
            Es.insertDataset(sdata, d)

        print Es
        #now fit this to MM equation to extract a Vmax for each [I]
        Es.fitDatasets(models=['Michaelis-Menten'], noiter=500, silent=True)
        for d in E.datasets:
            ferrs = Es.estimateExpUncertainty(d, runs=30)
            Es.addMeta(d, 'exp_errors', ferrs)

        v=[];iconcs=[]
        xerrors=[];yerrors=[]
        for d in Es.datasets:
            x = d.split('_')
            iconc = round(float(x[0]),3)
            vmax = Es.getMetaData(d)['Vmax']
            yerror = Es.getMeta(d, 'exp_errors')['Vmax'][1]
            v.append(vmax)
            iconcs.append(iconc)
            xerrors.append(iconc*self.sconcerr/100)
            yerrors.append(yerror)

        edata = EkinConvert.xy2ekin([iconcs, v], errorlists=(xerrors, yerrors),
                                        labels=('[I]','Vmax'))

        Es.insertDataset(edata, 'enzyme_conc')
        enzymeconc = self.getEnzymeConcfromAST(Es, dilution)
        #Es.saveProject('ast2test')
        return Es, enzymeconc


    def getEnzymeConcfromAST(self, E, dilution):
        """Extract enzyme conc from ast ekin project. Here we fit the
         velocities (or vmax) from either AST method to a line"""

        print 'finding enzyme conc. and error values..'
        d='enzyme_conc'
        edata=E.getDataset(d)
        ekd = EkinDataset(edata)
        E.fitDataset(d, model='Linear', conv=1e-12, silent=True)

        currfit = E.getFitData(d)
        expdata = Fitting.getExpData(edata)
        X = Fitting.makeFitter(currfit, expdata=expdata)
        currsqdiff = X.getpercSqDiff()
        try:
            newdata, newfit, X = self.fitportion(edata, currfit, percs=[0.8, 0.5])
        except:
            pass
        newsqdiff = X.getpercSqDiff()
        if newsqdiff < currsqdiff:
            E.data[d] = newdata
            E.setFitData(d, newfit)
            currsqdiff = newsqdiff

        ferr = E.estimateExpUncertainty(d, runs=50)
        E.addMeta(d, 'exp_errors', ferr)
        #print E.getMetaData(d)
        #Get fit and return [I] axis intercept
        a = E.getMetaData(d)['a']
        b = E.getMetaData(d)['b']
        I = - b/a
        #error is the sum square of relative errors in a and b
        erra = E.getMeta(d, 'exp_errors')['a'][1] / a
        errb = E.getMeta(d, 'exp_errors')['b'][1] / b
        #print a,b,erra,errb
        error =  abs(I * ( math.sqrt(pow(erra,2) + pow(errb,2)) ))
        print a,b,I, error, dilution
        enzymeconc = I * dilution
        error = error * dilution
        enzymeconc = (enzymeconc, error)
        return enzymeconc


    def checkFit(self, E, d, cutoff):
        """Rough check of any fit using percsqdiff and return status based on cutoff"""

        edata = E.data[d]
        ekd = EkinDataset(edata)
        currfit = E.getFitData(d)

        if currfit == None or len(currfit) == 0:
            return 0
        expdata = Fitting.getExpData(edata)
        X = Fitting.makeFitter(currfit, expdata=expdata)
        sqdiff = X.getpercSqDiff()
        #print d, sqdiff
        if sqdiff >= cutoff:
            return 0
        else:
            return 1

    def checkLinearFit(self, E, d, cutoff, fromstart=True, percs=None):
        """Check linear fit for raw abs data to get best fit.
            1. get best fit from 0.7, 0.5, 0.3 removed from end
            2. compare with original and pick the best fit
           We then return a flag stating if this point is to be used in the
           using percent sumsqdiff vals"""

        edata = E.data[d]
        ekd = EkinDataset(edata)
        currfit = E.getFitData(d)
        if currfit == None or len(currfit) == 0:
            return 0
        expdata = Fitting.getExpData(edata)
        X = Fitting.makeFitter(currfit, expdata=expdata)
        try:
            currsqdiff = X.getpercSqDiff()
        except:
            return 0
        currslope = X.getFitDict()['a']

        #remove different successive perc from end
        if fromstart==True:
            '''remove different successive perc from start'''
            newdata, newfit, X = self.fitportion(edata, currfit, percs=percs, X=X)
            newslope = X.getFitDict()['a']
            newsqdiff = X.getpercSqDiff()
            if newsqdiff <= currsqdiff and newslope > 0:
                E.data[d] = newdata
                E.setFitData(d, newfit)
                currsqdiff = newsqdiff
        else:
            newdata, newfit, X = self.fitportion(edata, currfit, percs=percs, X=X)
            newslope = X.getFitDict()['a']
            newsqdiff = X.getpercSqDiff()
            if newsqdiff <= currsqdiff and newslope > 0:
                E.data[d] = newdata
                E.setFitData(d, newfit)
                currsqdiff = newsqdiff

        #print 'perc fit error on point %s: %s' % (d, currsqdiff)
        if currsqdiff > cutoff:
            print 'sumsqdiff %s too high, exclude?' % currsqdiff
            return 0
        else:
            return 1


    def fitportion(self, edata, currfit, reverse=True, percs=None, X=None):
        """Compare fits using different portions of the data and return best.
           This method is not that reliable.. as its biased towards the fits with
           least no of points as they have lower sum sq diffs... """
        import copy
        if percs==None:
            percs = [0.7, 0.5, 0.3]

        currsqdiff = 100
        currdata = edata
        currX=X
        for perc in percs:
            tmpdata = copy.deepcopy(edata)
            ekd = EkinDataset(tmpdata)
            kys=ekd.DPkeys()
            lenp = int(len(kys)*perc)
            if reverse == True:
                i=len(kys)
                kys.reverse()
                for e in kys:
                    ekd.setDPActive(e, 0)
                    if i==lenp or not e in kys:
                        break
                    i-=1
            else:
                i=0
                for e in kys:
                    ekd.setDPActive(e, 0)
                    if i==lenp or not e in kys:
                        break
                    i+=1
            try:
                newfit, X = Fitting.doFit(tmpdata, model='Linear', conv=1e-8, grad=1e-8, silent=True)
                newsqdiff = X.getpercSqDiff()
            except:
                continue

            if newsqdiff < currsqdiff:
                best=perc
                currdata = tmpdata
                currfit = newfit
                currsqdiff = newsqdiff
                currX = X
        return currdata, currfit, currX


    def makeASTHtml(self):
        """make html page of ast raw fits"""
        ew = EkinWeb()
        path = os.path.join(self.path, 'plots')
        images=[]
        rownames=[]
        for name in self.required:
            E = EkinProject()
            E.openProject(os.path.join(self.path, name.replace(' ','')+'_ast'))
            E.datasets.sort()
            colnames=E.datasets
            rownames.append(name)
            for d in E.datasets:
                fname = os.path.join(path, name.replace(' ','')+d+'.png')
                E.plotDatasets(datasets=d, filename=fname, showfitvars=True, dpi=60)
                images.append(os.path.basename(fname))
        ew.makeHtml(images, outfile=os.path.join(self.path, 'ast.html'),
                            imgpath='plots', columns=9, title='AST results for set'+self.setname,
                            colheader=colnames, rowheader=rownames)
        return

    def makeKineticsHtml(self, run=None):
        """make html pages of kinetics raw and MM fits"""
        if run != None:
            runs=[run]
        else:
            runs=self.runs
        import gc
        ew = EkinWeb()
        path = os.path.join(self.path, 'plots')
        import copy
        rownames = copy.deepcopy(self.required)
        colnames = self.phlist
        sep='_'
        for r in runs:
            print 'doing run %s' %r
            kimages=[]
            mmimages=[]
            for name in self.required:
                gc.collect()
                E = EkinProject()
                E.openProject(os.path.join(self.path, name.replace(' ','')+'_'+str(r)))
                for ph in self.phlist:
                    dlist=[]
                    for s in self.subconcs:
                        dlist.append(ph+'_'+str(s))
                    fname = os.path.join(path, name.replace(' ','')+sep+str(r)+sep+ph+'_raw.png')
                    try:
                        E.plotDatasets(datasets=dlist, filename=fname, plotoption=3, legend=True, dpi=60)
                        kimages.append(os.path.basename(fname))
                    except:
                        kimages.append('dummy.png')

                    fname = os.path.join(path, name.replace(' ','')+sep+str(r)+'_mm_'+ph+'.png')
                    E.plotDatasets(datasets=ph, filename=fname, showfitvars=True, dpi=60)
                    mmimages.append(os.path.basename(fname))

            t='set '+self.setname+', run '+str(r)+' '+' raw kinetics fits'
            ew.makeHtml(kimages, outfile=os.path.join(self.path, 'run_'+str(r)+'_rawkinetics.html'),
                                imgpath='plots', columns=len(colnames), title=t,
                                colheader=colnames, rowheader=rownames)
            t='set '+self.setname+', run '+str(r)+' '+' mm fits'
            ew.makeHtml(mmimages, outfile=os.path.join(self.path, 'run_'+str(r)+'_MM.html'),
                                imgpath='plots', columns=len(colnames), title=t,
                                colheader=colnames, rowheader=rownames)
        return

    def makekcatsHtml(self, E=None, dataname='kcats'):
        """make html page of pka vs. kcat plots per variant"""
        if E==None:
            E = EkinProject()
            E.openProject(os.path.join(self.path, dataname))
        ew = EkinWeb()
        images=[]
        required = self.required
        runs = self.runs
        path = os.path.join(self.path, 'plots')
        if not os.path.exists(path):
            os.mkdir(path)
        colnames=['run 1','run 2','run 3','all runs']
        for name in required:
            tmp=[]
            for r in runs:
                fname = os.path.join(path, name.replace(' ','')+str(r)+dataname+'.png')
                E.plotDatasets(datasets=name+'_'+str(r), filename=fname, showfitvars=True, dpi=50)
                images.append(os.path.basename(fname))
                tmp.append(name+'_'+str(r))
            fname = os.path.join(path, name.replace(' ','')+dataname+'.png')
            E.plotDatasets(datasets=tmp, filename=fname, plotoption=3, legend=True,
                                 dpi=60)
            images.append(os.path.basename(fname))
        ew.makeHtml(images, outfile=os.path.join(self.path, dataname+'.html'),
                            imgpath='plots', columns=4, title='set '+self.setname,
                            colheader=colnames, rowheader=required)
        return

    def perReplicateHtml(self, data, errs, name, variants, setindex, savepath, runs):
        saveout = sys.stdout
        fsock = open(os.path.join(savepath,name+'_byreplicate.html'), 'w')
        sys.stdout = fsock
        print '<a>%s per replicate</a>' %name
        print '<table align=center width=60%>'
        print '<th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th>'\
            %('set', 'name','average','std dev.','error','rep 1','err','rep 2','err','rep 3','err')
        print '<tr>'
        for n in range(len(variants)):
            name = variants[n]
            kc=[]; kcerrs=[]
            if setindex != None:
                s = setindex[variants[n]]
            else:
                s = self.setname
            for r in runs:
                kc.append(data[r][n])
                kcerrs.append(errs[r][n])
            kcerr = numpy.sqrt(sum(numpy.square(kcerrs)));
            print '<td>%s</td> <td>%s</td> <td>%.2f</td> <td>%.2f</td> <td>%.2f</td>'\
                %(s, name, numpy.mean(kc), numpy.std(kc), kcerr)
            for r in runs:
                print '<td>%.2f</td> <td>%.2f</td>' %(data[r][n], errs[r][n])
            print '<tr>'
        print '</table>'
        return

    #plot summaries for variants
    def Summary(self, Ekc=None, Ekmkcat=None, Ekm=None, Ekslopes=None, path=None,
                        enzymeconcs=None, title=None, setindex=None, pathindex=None,
                        variants=None, width=8, height=18, axeslimits=False):
        import matplotlib
        import matplotlib.pyplot as plt
        import csv
        plt.rc("font", family="serif")
        #plt.rc('text', usetex=True)

        if Ekc==None and Ekmkcat==None:
            Ekc = EkinProject(); Ekmkcat=EkinProject(); Ekm=EkinProject();
            Ekslopes = EkinProject()
            Ekc.openProject(os.path.join(self.path, 'kcats'))
            Ekm.openProject(os.path.join(self.path, 'kms'))
            Ekmkcat.openProject(os.path.join(self.path, 'kcats_kms'))
            Ekslopes.openProject(os.path.join(self.path, 'kslopes'))
        totalvariants=0
        wtkcat=None;averagekcats=[]
        pkasfit=0
        dpkas={}; wtpka=None
        lines=[]

        if variants==None:
            variants = []
            for e in Ekc.datasets:
                #v=e.split('_')[0]
                v=e[:len(e)-2]
                if v not in variants:
                    variants.append(v)
        if enzymeconcs == None:
            enzymeconcs = self.enzymeconcs
        print 'variants', variants
        #runs = self.runs
        runs = [1,2,3]
        fig = plt.figure(figsize=(width, height))
        axesnames=['pkas','kcats','kms','kcat/km','deltakcat','kcat/wt']
        ax1 = fig.add_subplot(611, label=axesnames[0])
        ax2 = fig.add_subplot(612, label=axesnames[1])
        ax3 = fig.add_subplot(613, label=axesnames[2])
        ax4 = fig.add_subplot(614, label=axesnames[3])
        ax5 = fig.add_subplot(615, label=axesnames[4])
        ax6 = fig.add_subplot(616, label=axesnames[5])

        #get chosen wt data first
        #we use whatever wt is best here, or an average??
        if 'wt 1a_1' in Ekc.datasets:
            wtkcats=[]
            wtkms=[]
            wtpkas=[]; wtpkaerrs=[]
            for r in runs:
                d = 'wt 1a_'+str(r)
                ekd = EkinDataset(Ekc.getDataset(d))
                wtkcats.append(ekd.getYVal(7))
                ekmd = EkinDataset(Ekm.getDataset(d))
                wtkms.append(ekmd.getYVal(7))
                #if Ekc.getMeta(d, 'exp_errors').has_key('pKa'):
                    #wtpkas.append(Ekc.getMetaData(d)['pKa'])
                    #wtpkaerrs.append(Ekc.getMeta(d, 'exp_errors')['pKa'][1])
            wtkcat = numpy.mean(wtkcats)
            wtpka = numpy.mean(wtpkas)
            wtpkaerr = numpy.sqrt(sum(numpy.square(wtpkaerrs)))
            print self.bold+'wt kcat: '+self.reset, wtkcat
            print self.bold+'wt km: '+self.reset, numpy.mean(wtkms)


        kcats={};kcaterrs={};deltakcats={}
        deltakcaterrs={};pkas={}; pkaerrs={}
        kslopesall={}; kslopes={}; kcatskms={}
        kms={}; kmerrs={}; kmsall={}
        kcatskmsall={}

        print 'reading ekin projects'
        for r in runs:
            dpkas[r]={}
            pkas[r]=[]
            kcats[r]=[]
            kms[r]=[]
            kmsall[r]=[]
            kcatskms[r]=[]
            kcatskmsall[r]=[]
            kslopesall[r]=[]
            kslopes[r]=[]
            pkaerrs[r]=[]
            kcaterrs[r]=[]
            kmerrs[r]=[]
            kcatskmserrs=[]
            deltakcats[r]=[]
            deltakcaterrs[r]=[]
            kcatwtratios=[]
            enzconcs=[]
            lists = [kcats[r],kms[r],kcatskms[r],kcaterrs[r],kmerrs[r],kcatskmserrs,deltakcats[r],
                        deltakcaterrs[r],kcatwtratios,pkas[r],pkaerrs[r]]

            for name in variants:

                if path == None:
                    loadpath = self.path
                else:
                    #for combined we need to get the path for that variants' set
                    #print 'get path'
                    #print setindex[name]
                    loadpath = pathindex[setindex[name]]

                d = name+'_'+str(r)
                ed = Ekc.getDataset(d)
                if ed == None:
                    for l in lists:
                        l.append(-100)
                    continue
                totalvariants+=1
                ekd = EkinDataset(ed)

                if Ekc.getMeta(d, 'exp_errors').has_key('pKa'):
                    pKa = Ekc.getMetaData(d)['pKa']
                    pkas[r].append(pKa)
                    pkaerrs[r].append(Ekc.getMeta(d, 'exp_errors')['pKa'][1])
                    pkasfit+=1
                    if wtpka:
                        dpkas[r][name] = wtpka - pKa
                else:
                    pkas[r].append(0)
                    pkaerrs[r].append(0)
                    dpkas[r][name] = 0

                #we use ph 10 values for all results
                #delta kcat is now log of kcat mutant/wt
                kcat=ekd.getYVal(7)
                kcats[r].append(kcat)
                if wtkcat==None:
                    wtkcat=18
                xx = math.log(ekd.getYVal(7)/wtkcat)
                deltakcats[r].append(xx)
                deltakcaterrs[r].append(0.001)
                kcaterrs[r].append(ekd.getYError(7))
                ekd1 = EkinDataset(Ekm.getDataset(d))
                kms[r].append(ekd1.getYVal(7))
                kmerrs[r].append(ekd1.getYError(7))
                ekd2 = EkinDataset(Ekmkcat.getDataset(d))
                kcatskms[r].append(ekd2.getYVal(7))
                ekd3 = EkinDataset(Ekslopes.getDataset(d))
                kslopes[r].append(ekd3.getYVal(7))

                for val in ekd2.getY():
                    kcatskmsall[r].append(val)
                for val in ekd1.getY():
                    kmsall[r].append(val)
                kcatskmserrs.append(ekd2.getYError(7))
                if wtkcat:
                    kcatwtratios.append(kcat/wtkcat)
                enzconcs.append(Ekc.getMeta(d, 'enzymeconc'))

            ind=numpy.arange(len(variants))

            line1, = ax1.plot(ind, pkas[r], 'o', alpha=0.7)
            lines.append(line1)
            ax1.errorbar(ind, pkas[r], yerr=pkaerrs[r], fmt=None, ecolor=line1.get_color())
            line2, = ax2.plot(ind, kcats[r], 'o', alpha=0.7)
            ax2.errorbar(ind, kcats[r], yerr=kcaterrs[r], fmt=None, ecolor=line2.get_color())
            line3, = ax3.plot(ind, kms[r], 'o', alpha=0.7)
            ax3.errorbar(ind, kms[r], yerr=kmerrs[r], fmt=None, ecolor=line2.get_color())
            line4, = ax4.plot(ind, kcatskms[r], 'o', alpha=0.7)
            ax4.errorbar(ind, kcatskms[r], yerr=kcatskmserrs, fmt=None, ecolor=line2.get_color())
            line5, = ax5.plot(ind, deltakcats[r], 'o', alpha=0.7)
            ax5.errorbar(ind, deltakcats[r], yerr=deltakcaterrs[r], fmt=None, ecolor=line2.get_color())
            if wtkcat:
                line6, = ax6.plot(ind, kcatwtratios, 'o', alpha=0.7)

        #print 'kslopes', kslopes
        #print 'kcatskms', kcatskms

        print 'doing plots'
        axes = fig.get_axes()
        ft = matplotlib.font_manager.FontProperties(size=8)
        axes[0].legend(lines, runs, numpoints=1,loc='best')
        for a in axes:
            a.set_xticks(ind)
            a.set_xticklabels(variants, rotation='vertical',fontproperties=ft)
            a.grid(True)
            a.set_axis_bgcolor('#F5F5F5')

        ax1.set_ylim((4, 10))
        ax1.set_ylabel('pH')
        ax2.set_ylabel('$s^{-1}$')
        ax3.set_ylabel('[S](mM)')
        ax4.set_ylabel('$s^{-1}/mM$')
        ax5.set_ylabel('$\Delta kcat(s^{-1})$')
        ax5.set_ylabel('kcat/wt kcat')
        #to remove outliers that throw scaling off
        if axeslimits == True:
            ax2.set_ylim((0, 250))
            ax3.set_ylim((0, 5))
            ax4.set_ylim((0, 1500))
            ax5.set_ylim((0, 150))
            ax6.set_ylim((0, 5))

        ax1.set_title('fitted pka values per variant')
        ax2.set_title('kcat values @ ph 10')
        ax3.set_title('km values @ ph 10')
        ax4.set_title('kcat/km values @ ph 10')
        ax5.set_title('delta kcat values, ph 6-10')
        if wtkcat:
            ax6.set_title('kcat/wt kcat ratios')

        if path == None:
            savepath=self.path
        else:
            savepath=path
        fig.subplots_adjust(hspace=0.6)
        if title==None:
            title=self.setname+' plot summary'
        fig.suptitle(title, fontsize=20)
        fig.savefig(os.path.join(savepath, 'results.png'), dpi=150)

        #make html table and csv file of general stats
        colnames  = ['set', 'name', 'clonename', 'mutation', 'kcat', 'kcat err',
                     'dkcat' ,'dkcat err', 'km', 'kcat/km', 'kcat/km(slope)',
                      '[E]']
        cwriter = csv.writer(open(os.path.join(savepath,'results.csv'),'w'))
        cwriter.writerow(colnames)

        clonenames = self.getCloneNames()
        #print clonenames

        print 'found %s pkas out of %s clones/replicates' %(pkasfit, totalvariants)
        print 'writing stats'
        saveout = sys.stdout
        fsock = open(os.path.join(savepath,'stats.html'), 'w')
        sys.stdout = fsock
        self.doHeader()
        print '<table align=center width=60%>'
        for c in colnames:
            print '<th>%s</th>' %c
        print '<tr>'

        avgkcats=[]; avgkcaterrs=[]
        avgpkas=[]; avgpkaerrs=[]

        for n in range(len(variants)):
            name = variants[n]
            if not enzymeconcs.has_key(name):
                continue
            try:
                clonename = clonenames[name]['name']
                mut = clonenames[name]['mutation']
            except:
                clonename = ''
                mut = ''

            if type(enzymeconcs[name]) is tuple:
                e = enzymeconcs[variants[n]][0]
            else:
                e = enzymeconcs[variants[n]]
            if setindex != None:
                s = setindex[variants[n]]
            else:
                s = self.setname
            kc=[]; kcerrs=[]; dkc=[]; dkcerrs=[];
            km=[]
            kcatkm=[]; kslp=[]

            for r in runs:
                kc.append(kcats[r][n])
                kcerrs.append(kcaterrs[r][n])
                dkc.append(deltakcats[r][n])
                dkcerrs.append(deltakcaterrs[r][n])
                km.append(kms[r][n])
                kcatkm.append(kcatskms[r][n])
                kslp.append(kslopes[r][n])

            avgkcats.append(numpy.mean(kc))
            kcerr = numpy.sqrt(sum(numpy.square(kcerrs))); avgkcaterrs.append(kcerr)
            dkcerr = numpy.sqrt(sum(numpy.square(dkcerrs)))

            print '<td>%s</td> <td>%s</td> <td>%s</td> <td>%s</td> <td>%.2f</td> <td>%.2f</td> <td>%.2f</td> <td>%.2f</td> <td>%.2f</td> <td>%s</td> <td>%s</td>'\
                %(s, name, clonename, mut, numpy.mean(kc), kcerr, numpy.mean(dkc), dkcerr, numpy.mean(kcatkm), numpy.mean(kslp), e)
            print '<tr>'
            #write to csv also
            cwriter.writerow([s, name, clonename, mut, numpy.mean(kc), kcerr, numpy.mean(dkc), dkcerr, numpy.mean(km), numpy.mean(kcatkm), numpy.mean(kslp), e])
        print '</table>'
        sys.stdout.flush()
        sys.stdout = saveout
        fsock.close()

        #make html table of kcat per replicate stats
        self.perReplicateHtml(kcats, kcaterrs, 'kcats', variants, setindex, savepath, runs)
        self.perReplicateHtml(deltakcats, deltakcaterrs, 'kcats', variants, setindex, savepath, runs)
        self.perReplicateHtml(kms, kmerrs, 'kcats', variants, setindex, savepath, runs)

        sys.stdout.flush()
        sys.stdout = saveout
        fsock.close()
        #Plot of average kcats for all
        fig3 = plt.figure(figsize=(width,6))
        ax = fig3.add_subplot(111)
        w=1.0
        ind=numpy.arange(len(avgkcats))
        ax.bar(ind, avgkcats, width=w, yerr=avgkcaterrs, color='g', ecolor='r', lw=0.6,alpha=0.8)
        ax.set_xticks(ind+w/2)
        ax.set_xticklabels(variants, rotation='vertical',fontproperties=ft)
        fig3.suptitle('Average kcats', fontsize=20)
        ax.set_ylabel('$kcat (s^{-1})$')
        fig3.savefig(os.path.join(savepath, 'avgkcats.png'), dpi=150)

        #Bar plot of pkas for variants vs avg wt pka
        '''
        if wtpka!=None:
            fig4 = plt.figure(figsize=(width,6))
            ax = fig4.add_subplot(111)
            ind=numpy.arange(len(avgpkas))
            dpkas=[]; dpkaerrs=[]
            for i in range(len(avgpkas)):
                dpkas.append(avgpkas[i] - wtpka)
                dpkaerrs.append(avgpkaerrs[i] + wtpkaerr)
            ax.bar(ind, dpkas, width=w, yerr=dpkaerrs, color='b', ecolor='r', lw=0.6,alpha=0.8)
            ax.set_xticks(ind+w/2)
            ax.set_xticklabels(variants, rotation='vertical',fontproperties=ft)
            fig4.suptitle('$\Delta$pKa: wt - variant')
            ax.set_ylabel('$\Delta$pka')
            fig4.savefig('dpkas.png',dpi=150)'''

        return

    def compareKcatKmSlopes(self):
        """Compare kcatkm vs slopes"""
        print 'reading ekin projects'
        kcatskms={};kslopes={};kcatskmserrs={}
        for r in runs:
            kcatskms[r]=[]
            kslopes[r]=[]
            kcatskmserrs[r]=[]
            for name in variants:
                if path == None:
                    loadpath = self.path
                else:
                    loadpath = pathindex[setindex[name]]
                d = name+'_'+str(r)
                #get slope - kcat/km @ph10 from MM data
                try:
                    Emm = EkinProject()
                    Emm.openProject(os.path.join(loadpath, name.replace(' ','')+'_'+str(r)))
                    ekmm = EkinDataset(Emm.getDataset('10'))
                    kslopesph10[r].append(ekmm.getYVal(4))
                    #get all phs
                    for ph in self.phlist:
                        ekmm = EkinDataset(Emm.getDataset(ph))
                        kslopes[r].append(ekmm.getYVal(4))
                except:
                    kslopes[r].append(0)

        #Plot them
        fig2 = plt.figure(figsize=(10,14))
        ax1 = fig2.add_subplot(311)
        ax2 = fig2.add_subplot(312)
        ax3 = fig2.add_subplot(313)
        ll=[]
        for r in runs:
            diff=[]
            for k in range(len(kslopesall[r])):
                diff.append(abs(kslopesall[r][k]-kcatskmsall[r][k]))
            if len(kslopesph10[r]) != len(kcatskms[r]):
                continue
            a=ax1.plot(kslopesph10[r], kcatskms[r],'o', markersize=2.5)
            print len(kslopesall[r]), len(kcatskmsall[r])
            ax2.plot(kslopesall[r], kcatskmsall[r],'o', markersize=2.5)
            #ax3.plot(kmsall[r], diff, 'o', markersize=2.5)
            ll.append(a)
        cl = numpy.arange(0, max(kcatskmsall[r]))
        #ax2.plot(cl, cl, 'r', alpha=0.5)
        ax1.set_xlabel('slope')
        ax1.set_ylabel('fitted kcat/km')
        ax1.legend(ll, runs, numpoints=1,loc='best')
        ax1.set_title('@ph10'); ax2.set_title('all ph')
        ax3.set_title('km vs diff')
        fig2.suptitle('slope vs kcat/km from fit', fontsize=20)
        fig2.savefig(os.path.join(savepath,'kcatskmsvsslope.png'), dpi=150)

        return

    def comparedpKas(self):
        """predicted dpka vs. our crappy dpkas"""
        import matplotlib.pyplot as plt
        f= open('/enc/Novozymes/Data/Hwbb/jan/dpkas.csv')
        cw = csv.reader(f)
        predicted={}
        for r in cw:
            name = string.lower(r[0]+' '+r[1])
            predicted[name]=r[2]

        fig=plt.figure()
        ax=fig.add_subplot(111)
        for r in self.runs:
            xd=[];yd=[]
            for p in predicted:
                if p in self.dpkas[r]:
                    xd.append(self.dpkas[r][p])
                    yd.append(predicted[p])
            ax.plot(xd,yd,'x')
        ax.set_title('dpKas: fitted vs. predicted')
        ax.set_xlabel('fitted experimental')
        ax.set_ylabel('predicted')
        fig.savefig('dpkasvspredicted.png',dpi=100)
        return

    def readEnzymeConcs(self):
        """Read previously calculated enzyme concs from the file"""
        ef = open(os.path.join(self.path,'enzymeconcs.txt'))
        cw=csv.reader(ef)
        enzconcs = {}
        templateconcs = {}
        for r in cw:
            enzconcs[r[0]] = float(r[1])
            templateconcs[r[0]] =  float(r[3])
        #print enzconcs
        return enzconcs, templateconcs

    def enzymeConcsfromTemplate(self):
        """Use enz concs from template"""
        print 'Using enzyme concentrations in %s' %self.astdilfile
        print
        self.astdilutions, xx = self.getDilutions(os.path.join(self.path, self.astdilfile))
        enzconcs={}
        for n in self.template_enzconcs:
            enzconcs[n] = (self.template_enzconcs[n], 0.01)
        print 'template ast concentrations: ', enzconcs
        return enzconcs

    def temperatureResults(self, ph='7', sets=None):
        """Temp data"""

        symb = unichr(177).encode("UTF-8")
        import csv
        clonenames = self.getCloneNames(allsets=False)
        dups=['1','2']
        tms=[]; labels=[]

        #get wt Tm first.. pick best ones
        wtlist = ['wt 1a','wt 2a','wt 4a','wt 2','wt 3','wt 5']
        Ewt = EkinProject()
        pth1 = 'Novozymes/Data/Hwbb/march/setJ/'
        pth2 = 'Novozymes/Data/Hwbb/jan/setF/'
        wtvals = []
        for name in wtlist:
            print name
            if name in ['wt 1a','wt 2a','wt 4a']:
                Ewt.openProject(os.path.join(pth1, name.replace(' ','')+'_temp'))
            else:
                Ewt.openProject(os.path.join(pth2, name.replace(' ','')+'_temp'))
            wtvals.append(Ewt.getMetaData('1_'+ph)['t50'])
            #wtTmerr = Ewt.getMeta(d, 'exp_errors')['t50'][1]

        wtTm = numpy.mean(wtvals)
        wtTmerr = round(numpy.std(wtvals),2)
        print 'wt Tm:', wtTm, symb, wtTmerr

        def plotDups(x, labels):
            a=[i[0] for i in x]
            b=[i[1] for i in x]
            import matplotlib.pyplot as plt
            #plt.rc('text', usetex=True)
            fig=plt.figure()
            ax=fig.add_subplot(111)
            ax.plot(a,b,'x',alpha=0.7)
            ax.plot(a,a,'-')
            ax.set_title('Tm consistency between plate duplicates')
            ax.set_xlabel('Tm1')
            ax.set_ylabel('Tm2')
            for i in x:
                if abs(i[0]-i[1])>=5:
                    j=x.index(i)
                    ax.text(i[0],i[1]+0.5,labels[j],size=8)
            fig.savefig('dupTm_correl.png',dpi=100)
            return

        def writeresults(outpath=self.path):
            for name in self.required:
                E=EkinProject()
                E.openProject(os.path.join(self.path, name.replace(' ','')+'_temp'))
                t=[]; dt=[]
                i=0
                for dup in dups:
                    d = dup+'_'+ph
                    fdata = E.getFitData(d)
                    if not d in E.datasets: continue
                    tm = E.getMetaData(d)['t50']
                    tmerror = E.getMeta(d, 'exp_errors')['t50'][1]
                    dt.append(round(wtTm-tm,3))
                    try:
                        clonename = clonenames[name]['name']
                        mut = clonenames[name]['mutation']
                    except:
                        clonename = ''; mut = ''
                    t.append(round(tm,2))
                if not d in E.datasets: continue

                avdeltatm = round(numpy.mean(dt),2)
                cw.writerow([self.setname, name, clonename, mut, t[0], t[1], dt[0], dt[1], avdeltatm, wtTmerr])
                print name, clonename, mut, t[0], t[1], dt[0], dt[1], avdeltatm, wtTmerr
                tms.append(t)
                labels.append(name)
            return
        p = ph.replace(',','.')
        header = ['set', 'name', 'clonename', 'mut', 'tm1_'+p, 'tm2_'+p, 'dtm1_'+p, 'dtm2_'+p, 'avdeltatm_'+p,'error']
        if sets == None:
            cw = csv.writer(open(os.path.join(self.path,'temperature.csv'),'w'))
            cw.writerow(header)
            writeresults()
        else:
            setfiles=[]
            cw = csv.writer(open(os.path.join('.','temperature_%s.csv' %p),'w'))
            cw.writerow(header)
            for s in sets:
                setfiles.append('set'+s+'.defs')
            for s in setfiles:
                self.parseDefs(s)
                writeresults()
        plotDups(tms, labels)
        #print tms
        return

    def doResults(self, Ekcats=None, all=False):
        """Make results summary/plots"""
        print 'producing results summaries..'
        if all == True:
            print 'doing kcat etc. hmtl plots..'
            self.makekcatsHtml()
            self.makekcatsHtml(dataname='kms')
            self.makekcatsHtml(dataname='kcats_kms')
            self.makekcatsHtml(dataname='kslopes')

        enz,templ=self.readEnzymeConcs()
        if len(enz) == 0:
            enz = self.enzymeConcsfromTemplate()
        self.Summary(Ekcats, enzymeconcs=enz)

        return

    def combinedSummary(self, sets=None, variants=None):
        """Put all results in one project and plot summary"""
        if sets == None:
            sets = self.sets
        setfiles=[]
        for s in sets:
            setfiles.append('set'+s+'.defs')
        if variants!=None:
            print ' doing variants', variants
        Ekcat=EkinProject()
        Ekcatkm=EkinProject()
        Ekm=EkinProject()
        Ekslopes=EkinProject()
        enzconcs={}
        templconcs = {}
        setindex = {}
        pathindex = {}
        conflict=[]
        for s in setfiles:
            self.parseDefs(s)
            E1=EkinProject()
            E2=EkinProject()
            E3=EkinProject()
            E4=EkinProject()
            E1.openProject(os.path.join(self.path, 'kcats'))
            E2.openProject(os.path.join(self.path, 'kcats_kms'))
            E3.openProject(os.path.join(self.path, 'kms'))
            E4.openProject(os.path.join(self.path, 'kslopes'))
            Ekcat.addProject(E1)
            Ekcatkm.addProject(E2)
            Ekm.addProject(E3)
            Ekslopes.addProject(E4)
            try:
                enz,templ=self.readEnzymeConcs()
                templconcs.update(templ)
            except:
                enz = self.enzymeConcsfromTemplate()

            pathindex[self.setname] = self.path
            for v in self.variantnames:
                #if not v in enzconcs.keys():
                setindex[v] = self.setname

            for k in enz:
                if not k in enzconcs.keys():
                    enzconcs[k] = enz[k]
                else:
                    pass

        if variants == None:
            w=20
        else:
            w=int(math.sqrt(len(variants)))+2
        self.Summary(Ekcat, Ekcatkm, Ekm, Ekslopes, path=os.getcwd(), enzymeconcs=enzconcs,
                            title='all sets plot summary',
                            width=w, height=18, axeslimits=True,
                            setindex=setindex, variants=variants,
                            pathindex=pathindex)
        return

    def plotTemplatevsenzConcs(self, sets):
        """Compare our enzyme concs vs template, testing purposes only"""
        import matplotlib
        import matplotlib.pyplot as plt
        setfiles=[]
        for s in sets:
            setfiles.append('set'+s+'.defs')
        enzconcs={}
        templconcs = {}
        pathindex = {}
        for s in setfiles:
            self.parseDefs(s)
            try:
                enz,templ=self.readEnzymeConcs()
                templconcs.update(templ)
            except:
                enz = self.enzymeConcsfromTemplate()
            pathindex[self.setname] = self.path
            for k in enz:
                if not k in enzconcs.keys():
                    enzconcs[k] = enz[k]
                else:
                    pass

        if enzconcs !=None:
            print 'plotting enzyme concs. ours vs template'
            fig1=plt.figure()
            ax=fig1.add_subplot(111)
            xd=[]; yd=[]
            for n in enzconcs:
                if enzconcs[n]<20 and enzconcs[n]>0:
                    xd.append(templconcs[n])
                    yd.append(enzconcs[n])
            ax.plot(xd, yd, 'x')
            ax.plot(xd, xd, '-')
            ax.set_title('enzyme concentrations from AST: PEAT vs template')
            ax.set_ylabel('PEAT')
            ax.set_xlabel('template')
            ax.grid(True)
            fig1.savefig('enzymecvstemplate.png')
        return


    def doAST(self, noraw=False):
        """Do AST using current defs
           Ekinraw is the optional manually re-fitted raw ast data"""
        enzymeconcs={}
        enzymeconcs1={}
        self.astdilutions, xx = self.getDilutions(os.path.join(self.path, self.astdilfile))
        if hasattr(self, 'template_enzconcs'):
            print self.template_enzconcs

        self.astraw = self.loadASTData(path=self.path, variantnames=self.variantnames, setname=self.setname,
                                        date=self.astdate, inhibconcs=self.subconcs)

        for name in self.required:
            print 'processing ast for variant', name

            #if noraw is true we don't get the raw data, but use the already made
            #ekin fits to extract the enzyme concs - handy for manually refit ast data
            if noraw == True:
                Ea=EkinProject()
                Ea.openProject(os.path.join(self.path, name.replace(' ','')+'_ast'))
                enzymeconc = self.getEnzymeConcfromAST(Ea, dilution=self.astdilutions[name])
            else:
                Ea, enzymeconc = self.fitAST(self.astraw[name], dilution=self.astdilutions[name])
                #Ea, enzymeconc1 = self.fitAST2(self.astraw[name], dilution=self.astdilutions[name])
                Ea.saveProject(os.path.join(self.path, name.replace(' ','')+'_ast'))

            enzymeconcs[name]=enzymeconc
            #enzymeconcs1[name]=enzymeconc1

        cw =  csv.writer(open(os.path.join(self.path, 'enzymeconcs.txt'),'w'))
        print 'found enzyme concentrations:'
        print '-----------------------------------'

        symb = unichr(177).encode("UTF-8")
        for name in enzymeconcs:
            print self.bold + name + self.reset, enzymeconcs[name][0], symb, enzymeconcs[name][1], '(%s)' %self.template_enzconcs[name]
            #print self.bold + name + self.reset, enzymeconcs1[name][0], symb, enzymeconcs1[name][1], '(%s)' %self.template_enzconcs[name]
            cw.writerow([name, enzymeconcs[name][0], enzymeconcs[name][1], self.template_enzconcs[name]])

        print '-----------------------------------'
        self.enzymeconcs = enzymeconcs

        E=EkinProject()
        y = [enzymeconcs[i][0] for i in enzymeconcs]
        x = [self.template_enzconcs[i] for i in enzymeconcs]
        yerrs = [enzymeconcs[i][1] for i in enzymeconcs]
        data = EkinConvert.xy2ekin([x, y], errorlists=(None, yerrs), labels=('templ','ours'))
        E.insertDataset(data, 'templ. vs simple', replace=True)
        #y = [enzymeconcs1[i][0] for i in enzymeconcs1]
        #data = EkinConvert.xy2ekin([x, y], labels=('templ','ours'))
        #E.insertDataset(data, 'templ. vs complex', replace=True)
        E.fitDatasets(models=['Linear'])
        E.saveProject(os.path.join(self.path, 'astvstemplate'))

        return enzymeconcs


    def doKinetics(self, enzymeconcs=None, noraw=False, doftest=False):
        """Do Kinetics using current settings"""
        if enzymeconcs == None:
            enzymeconcs=self.enzymeconcs
        required = self.required
        runs = self.runs
        ekinprojects={}
        self.kslopesph10={}

        for r in runs:
            print 'kinetics data, doing run %s' %r
            #remove old log file
            if os.path.exists(os.path.join(self.path,'log_'+str(r)+'.html')):
                os.remove(os.path.join(self.path,'log_'+str(r)+'.html'))

            #check if we have alt dilutions for this run, or use default ones
            if self.otherdilutions != '' and self.otherdilutions.has_key(r):
                f = self.otherdilutions[r]
                print 'Using alternate dilutions for run %s' %r
                rawdilutions, rawdilutions_alt  = self.getDilutions(os.path.join(self.path, f))
                #print rawdilutions
            else:
                rawdilutions, rawdilutions_alt = self.getDilutions(os.path.join(self.path, self.rawdilfile))

            if noraw == False:
                raw = self.loadRawData(path=self.path, variantnames=self.variantnames, phvalues=self.phlist, run=r,
                                      date=self.runs[r], setname=self.setname, subconcs=self.subconcs)

            for name in required:
                if name == 'empty':
                    continue
                print 'processing variant %s with enzyme conc. %s' %(name, enzymeconcs[name])
                if rawdilutions_alt.has_key(name):
                    altdilution = rawdilutions_alt[name]
                else:
                    altdilution = None
                if noraw == False:
                    Ek = self.processKineticsRaw(raw[name], phvalues=self.phlist)
                    Ek.saveProject(os.path.join(self.path, name.replace(' ','')+'_'+str(r)))
                Ek = self.analyseKinetics(name=name.replace(' ','')+'_'+str(r),
                                            enzymeconc=enzymeconcs[name][0],
                                            enzymeconcerr=enzymeconcs[name][1],
                                            dilution=rawdilutions[name],
                                            altdilution=altdilution,
                                            doftest=doftest,
                                            run=r)
                Ek.saveProject(os.path.join(self.path, name.replace(' ','')+'_'+str(r)))
                ekinprojects[name] = Ek

        print 'finished.. '
        return Ek

    def getCloneNames(self, allsets=True):
        """Get map of clone/seq names to plate/well name used here, we now use
           the master list stored in novo db labbook"""

        import PEATSA.Core as Core
        pdb = '/enc/10R_TS.pdb'
        confseq = self.DB.getLabbookData('confirmedseq',
                        ['plate','well','Tapo.check','variant'])
        seqnames={}
        for c in zip(confseq['plate'],confseq['well'],confseq['Tapo.check'],confseq['variant']):
            name = (c[0]+' '+c[1]).lower()
            if c[2] != '':
                #convert from tapo format..
                m = c[2]
                if '*' in m: continue
                x = m.split(' ')
                if '' in x: x.remove('')
                mut=''
                for i in x: mut+='A'+i+'+'
                mut=mut[:-1]
                s = Core.Data.MutationSet(code=mut)
                seqnames[name] = (mut, c[3])

        allvariants = []
        if allsets == True:
            for s in self.sets:
                self.parseDefs(os.path.join(self.defsfilepath, 'set'+s+'.defs'))
                allvariants.extend(self.variantnames)
        else:
            allvariants = self.variantnames
        names = {}
        cw = csv.writer(open('clonenames.csv','w'))

        for variant in seqnames:
            name = seqnames[variant][1].strip('NZo:')
            names[variant] = {'name': name}
            mutation = seqnames[variant][0]
            names[variant]['mutation'] = mutation
            if variant in allvariants:
                #print name, variant, mutation
                cw.writerow([name, variant, mutation])

        return names

    def sendtoPEAT(self, server='localhost'):
        """Send all the data to PEATDB"""
        from PEATDB.Base import PDatabase
        self.DB = PDatabase(server=server, username='farrell',
                         password='123', project='novo')
        print self.DB
        DB=self.DB
        DB.addField('setname','text')
        DB.addField('clonename','text')
        #DB.addField('mutation','text')
        DB.addField('ast','General')
        DB.addField('enz_conc','text')
        DB.addField('temp','General')
        DB.addField('t50','text')
        DB.addField('cd','General')

        clonenames = self.getCloneNames(allsets=False)
        enzymeconcs = self.enzymeConcsfromTemplate()
        kprj = {'kcats':EkinProject(), 'kcats_kms':EkinProject(), 'kslopes':EkinProject()}
        for k in kprj:
            kprj[k].openProject(os.path.join(self.path, k))

        for name in self.required:
            DB.add(name)
            DB[name].setname = self.setname
            print name
            if clonenames.has_key(name):
                DB[name].clonename = clonenames[name]['name']
                DB[name].Mutations = clonenames[name]['mutation']
                
            #kinetics data              
            DB[name]['enz_conc'] = str(enzymeconcs[name][0])
            Ea = EkinProject()
            Ea.openProject(os.path.join(self.path, name.replace(' ','')+'_ast'))
            DB[name]['ast'] = Ea
            
            for r in self.runs:
                print 'adding run %s data' %r
                E = EkinProject()
                E.openProject(os.path.join(self.path, name.replace(' ','')+'_'+str(r)))
                col = 'raw_'+str(r)
                DB.addField(col,'General')
                DB[name][col] = E

                #add kcat etc. vals for each run
                '''for k in kprj:
                    for r in self.runs:
                        d = name+'_'+str(r)
                        col = k+'_'+str(r)
                        ekd = EkinDataset(kprj[k].getDataset(d))
                        DB.addField(col,'text')
                        DB[name][col] = str(round(ekd.getYVal(7),5))'''
            

            #temp data
            Et = EkinProject()
            Et.openProject(os.path.join(self.path, name.replace(' ','')+'_temp'))
            try:
                t50=numpy.mean(Et.getMetaData('1_7')['t50'] + Et.getMetaData('2_7')['t50'])
                DB[name]['t50'] = round(t50,2)
                DB[name]['temp'] = Et
                DB.updateCellCache(name, 't50')
                DB.updateCellCache(name, 'temp')
            except:
                pass

        print DB        
        DB.commit(note='added from script')
        return DB

    def createOptions(self):
        from optparse import OptionParser
        parser = OptionParser()

        parser.add_option("-d", "--defs", dest="defs",
                            help="Provide a defs file with setup", metavar="FILE")
        parser.add_option("-e", "--ekinast", dest="ekinast", action='store_true',
                            help="Calculate enzyme concentrations from current ekin ast prjs", default=False)
        parser.add_option("-z", "--useenzymeconcs", dest="useenzymeconcs", action='store_true',
                            help="Use enzyme concentrations from a file", default=False)
        parser.add_option("-m", "--ekinmm", dest="ekinmm",  action='store_true',
                            help="Use current raw mm calculations from ekin projects", default=False)
        parser.add_option("-n", "--nocalcs", dest="nocalcs",  action='store_true',
                            help="Don't do any calculations", default=False)
        parser.add_option("-f", "--doftest", dest="doftest",  action='store_true',
                            help="Do f-test for pka fits", default=False)
        parser.add_option("-p", "--plots", dest="plots",  action='store_true',
                            help="Make html plots of kcat data ", default=False)
        parser.add_option("-w", "--rawplots", dest="rawplots",  action='store_true',
                            help="Make html plots of raw/ast data ", default=False)
        parser.add_option("-c", "--combine", dest="combine",  action='store_true',
                            help="Combine multiple set results into one plot", default=False)
        parser.add_option("-g", "--debug", dest="debug",  action='store_true',
                            help="Debug set true, prints out more info", default=False)
        parser.add_option("-v", "--variants", dest="variants", default=None,
                            help="Select variants for combined results", metavar="FILE")
        parser.add_option("-l", "--run", dest="run", default=None,
                            help="Specify a replicate")
        parser.add_option("-k", "--dokcatproj", dest="dokcatproj", action='store_true',
                            help="Skip to creating kcat/km projects from current", default=False)
        parser.add_option("-r", "--results", dest="doresults", action='store_true',
                            help="Just create results", default=False)
        parser.add_option("-t", "--test", dest="test", action='store_true',
                            help="Just test", default=False)
        parser.add_option("-x", "--clonenames", dest="clonenames", action="store", type="string",
                            help="Get clone names", default=None,)
        parser.add_option("-s", "--peat", dest="peat", action='store_true',
                            help="Send current set to peat", default=False)
        return parser

    def run(self, opts):
        """Application to test the class"""
        app = self
        if opts.defs != None:
            app.parseDefs(opts.defs)
            if opts.nocalcs == True:
                print 'defs file ok, not doing calcs'
                print
                if opts.plots == True:
                    app.doResults(all=False)
            elif opts.test == True:
                #enzymeconcs = app.doAST() #test stuff here
                #app.makeASTHtml()
                #app.plotTemplatevsenzConcs(sets=['B','C','D','E','F'])
                #raw = app.loadTempData()
                #app.processTempData(raw)
                app.temperatureResults(ph='7', sets=['B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'])

            elif opts.peat == True:
                DB = app.sendtoPEAT()
            elif opts.doresults == True:
                app.doResults()
            else:
                if opts.dokcatproj == True:
                    app.combineKcatProjects()
                else:
                    if opts.useenzymeconcs == True:
                        enzymeconcs = app.enzymeConcsfromTemplate()
                    elif opts.ekinast == True:
                        enzymeconcs = app.doAST(noraw=True)
                    else:
                        enzymeconcs = app.doAST(noraw=False)
                    app.doKinetics(enzymeconcs, noraw=opts.ekinmm, doftest=opts.doftest)
                    app.combineKcatProjects()
                app.doResults()

            if opts.rawplots == True:
                print 'producing raw html plots..'
                app.makeASTHtml()
                app.makeKineticsHtml(run=opts.run)

        elif opts.combine != False:
            print 'combining sets results'
            v=None
            if opts.variants != None:
                v=[]
                fv=open(opts.variants ,'r')
                for r in fv.readlines():
                    n=r.strip('\n')
                    if n!='':
                        v.append(n)
            app.combinedSummary(variants=v)
        elif opts.clonenames != None:
            app.getCloneNames()
        elif opts.peat == True:
            #for s in app.sets:
            for s in ['H','I','J','K','L','M','N','O','P']:
                app.parseDefs('set'+s+'.defs')
                DB = app.sendtoPEAT()
            #DB.createCellCache()
            #DB.commit()
        else:
            print 'No defs file found. Use -d and provide one'
        return

if __name__ == '__main__':
    try:
        from PEATDB.Base import PDatabase
        DB = PDatabase(server='localhost', username='farrell',
                         password='123', project='novo')
    except:
        DB=None
    print DB
    app = KineticsAnalyser(db=DB)
    opts, remainder = app.parser.parse_args()
    app.run(opts)

