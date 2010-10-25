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

import cgi
import cgitb
#cgitb.enable()
import math
import sys, os, types
import StringIO
import tempfile
from PEATDB.Ekin.Base import EkinProject

class EkinWeb():
    """Class for opening and displaying Ekin file datasets.
       This is a server side script for providing html tables of ekin plots
       and possibly some other html representations of ekin processed data
       Uses the EkinProject class to do the plots.
    """

    def __init__(self, DB=None):
        """ """
        self.DB = DB    #we pass the peat DB
        self.columns=2
        self.graphwidth=6
        self.graphheight=4
        self.fontsize=12
        self.labelsize=8
        self.pointsize=4.0
        self.legend=0
        #self.doheader()
        tempfile.tempdir=os.getcwd()

        return

    def doheader(self, title):
        print '<head>'
        print '<title>%s</title>' %title
        #print '<link href="http://peat.ucd.ie/titration_db/styles.css" rel="stylesheet" type="text/css" />'
        print '</head>'
        return

    def dofooter(self):
        """Footer"""
        print '<br>'
        return

    def maketempImage(self):
        imgfile = tempfile.mktemp('.png')
        return imgfile

    def makeHtml(self, images, outfile=None, imgpath=None, columns=2,
                        title=None,colheader=None, rowheader=None):
        """Put given images into a html file and save it"""
        if outfile != None:
            saveout = sys.stdout
            fsock = open(outfile, 'w')
            sys.stdout = fsock
        if title!=None:
            print '<h1>%s</h1>' %title
        print '<table id="mytable" align=center cellspacing="0" borderwidth=1>'
        if colheader!=None:
            print '<tr>'
            print '<th></th>'
            for col in range(columns):
                print '<th>%s</th>' %colheader[col]
            print '</tr>'
        row=0;c=0
        for img in images:
            if c==0:
                print '<tr>'
                if rowheader!=None:
                    try:
                        print '<th>%s</th>' %rowheader[row]
                    except:
                        print '<th>??</th>'

            print '<td> <img src=%s/%s  align=center></td>' % (imgpath,img)
            c=c+1
            if c >= columns:
                print '</tr>'
                row=row+1
                c=0
        print '</table>'
        if outfile != None:
            sys.stdout.flush()
            sys.stdout = saveout
            fsock.close()
        return

    def showEkinPlots(self, ekindata=None, project=None, filename=None,
                            datasets='ALL',
                            title='Ekin plots',
                            outfile=None, imgpath=None, path='',
                            normalise=False, showfitvars=False,
                            plotoption=1, columns=2, legend=False, size=(8,6),
                            logx=False, logy=False):
        """Plot ekin datasets from the provided ekin project data"""

        def downloadLink():
            #do csv download link for data displayed
            print '<table id="mytable" valign=top>'
            cfname = tempfile.mktemp('.csv', dir=csvpath)
            E.exportCSV(filename=cfname)
            p = os.path.split(path)[0]
            print '<td><a href=%s title="right-click to save as"> download data </a>' \
                      %os.path.join(p, 'csv', os.path.basename(cfname))
            print '</td></tr>'
            print '</table>'

        csvpath = os.path.join( os.path.split(imgpath)[0], 'csv')

	print project
        if ekindata != None: #convert from ekindata
            E = EkinProject(data=ekindata, mode='NMR titration')
        elif project != None: #just passed object
            E = project
        elif filename != None: #load project from file
            E = EkinProject()
            E.openProject(project)
        else:
            return
	E.checkDatasets()
        #if outfile is given, we override imgpath
        if outfile != None and imgpath==None:
            imgpath = os.path.dirname(outfile)
        if imgpath != None:
            tempfile.tempdir = imgpath

        size=(8,6)
        if datasets == 'ALL':
            #we plot all the datasets
            datasets = E.datasets
        if plotoption == 1:
            if columns>2:
                size=(4,3)
            imagenames={}
            for d in datasets:
                imgfile = self.maketempImage()
                name = os.path.basename(imgfile)                
                E.plotDatasets(d, filename=imgfile, size=size, linecolor='r',
                                    normalise=normalise,showfitvars=showfitvars,legend=legend,
                                    logx=logx, logy=logy)
                imagenames[d] = name

        elif plotoption == 3:
            name = self.maketempImage()
            E.plotDatasets(datasets, filename=name, plotoption=3, size=size,
                                normalise=normalise, legend=legend,
                                logx=logx, logy=logy)

        if outfile != None:
            saveout = sys.stdout
            fsock = open(outfile, 'w')
            sys.stdout = fsock

        self.doheader(title)
        downloadLink()

        print '<table id="mytable" align=center cellspacing="0" borderwidth=1>'
        row=1;c=1
        datasets.sort()
        if plotoption == 1:
            for d in datasets:
                if not imagenames.has_key(d):
                    continue
                if c==1:
                    print '<tr>'
                print '<td> <img src=%s/%s  align=center></td>' % (path, imagenames[d])
                print '<td class="alt">'
                #use ekinproject to supply formatted fit and model info here..
                self.showMetaData(ekinproj=E, dataset=d)
                print '</td>'
                c=c+1
                if c >= columns:
                    print '</tr>'
                    row=row+1
                    c=1
        elif plotoption == 3:
            print '<td> <img src=%s/%s  align=center></td>' % (path, os.path.basename(name))
            print '<td class="alt">'
            #use ekinproject to supply formatted fit and model info here..
	    x=1
            for d in datasets:
		if x>2:
		   n=True
		   x=0
	    	else:
		   n=False
		if n==False:
		   print '<td>'
            	self.showMetaData(ekinproj=E, dataset=d)
		x+=1
            print '</td>'
            c=c+1
            if c >= columns:
                print '</tr>'
                row=row+1
                c=1


        print '</table>'
        if outfile != None:
            sys.stdout.flush()
            sys.stdout = saveout
            fsock.close()
        return

    #@classmethod
    def showMetaData(self, ekinproj=None, ekindata=None, dataset=None, fdata=None, silent=False):
        """Print html of fit and metadata for the given dataset"""

        if fdata == None:
            if ekinproj == None and ekindata!=None:
                E = EkinProject(data=ekindata)
            else:
                E = ekinproj

            fdata = E.getMetaData(dataset)
        fsock = None
        if silent == True:
            saveout = sys.stdout
            sys.stdout = fsock = StringIO.StringIO()

        print '<table id="mytable">'
        print '<tr>'
        print '<td class="alt" style="bold" colspan=2>%s</td><tr>' % dataset
        kys = fdata.keys()
        ignore = ['error']
        for k in sorted(fdata):
            if k in ignore:
                continue
            if fdata[k] == None:
                continue
            elif type(fdata[k]) is types.DictType:
                print '<td>%s</td>' %k
                print '<td> <table id="mytable">'
                for n in fdata[k]:
                    val = fdata[k][n][1]
                    print '<td class="alt">%s</td><td>%.2f</td><tr>' %(n, val)
                print '</table></td><tr>'
            elif type(fdata[k]) is types.StringType:
                print '<td class="alt">%s</td><td>%s</td><tr>' %(k, fdata[k])

            else:
                print '<td class="alt">%s</td><td>%.2f</td><tr>' %(k, fdata[k])
        print '</table>'
        if silent == True:
            sys.stdout = saveout

        if fsock == None:
            return ''
        else:
            return fsock.getvalue()
        return

