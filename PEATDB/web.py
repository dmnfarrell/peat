#!/bin/env python
#
#
# This file is part Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-
# All rights reserved
#
#
# Rewritten by D Farrell, Jan 2009
#

import cgi
import sys,os, string, types
os.environ['MPLCONFIGDIR']='/tmp'
from PEATDB.Base import PDatabase
from PEATDB.Record import PEATRecord
from PEATDB.PEATTables import PEATTableModel, PEATTable

class PEATWeb:

    """Class providing a web interface for PEAT projects"""

    def __init__(self, server='localhost', project='test', port=8080,
                        user=None, passwd=None,
                        bindir='', fullpath=''):
    
        """bindir : path to cgi script in server address
           fullpath : file system path to script  """

        import socket
        self.host = socket.getfqdn(socket.gethostname())
      
        self.sessionkey=None
        self.form=cgi.FieldStorage()
        self.server = server
        self.project = project
        self.user = user
        self.password = passwd
        self.port = port
        
        self.bindir = bindir
        dirname = os.path.basename(fullpath)
        self.imgdir = '/' + dirname + '/images'
        self.plotsdir = '/'+ dirname + '/plots'
        self.imagepath = fullpath + '/plots'

        action = self.form.getfirst('action')
        self.action = action
        if not self.form.getfirst('action'):
            self.show_intro()
            return
        elif action == 'show_intro':
            self.show_intro()
        elif action == 'show_help':
            self.show_help()
        elif action == 'show_all':
            self.show_all()
        elif action=='show_specific':
            self.show_specific()
        elif action == 'show_datasets':
            self.show_datasets()
        elif action == 'search':
            self.show_search_results()

        return

    def show_intro(self):
        '''Show intro page'''
        self.show_DB_header(menu=1)        
        print '<div align=left>'
        print '<big><big><a>PEAT DB Web Interface, %s </big></a>' %self.project
        print '<br>'
        print '<UL>\
                <LI>Browse contents of the DB in one table\
                <LI>Visualise and overly any combination of curves (and our fits)\
                <LI>Download the raw data as csv/text files and refit/analyse\
                </UL>'
        print '</div>'
        self.footer()
        return

    def show_help(self):
        '''Show help page'''
        self.show_DB_header(menu=1)

        print '<div align=left>'
        print '<h2>Help information for using these pages is available <a href="http://enzyme.ucd.ie/main/index.php/"> here</a></h2>'
        print '</div>'
        self.footer()
        return

    def show_DB_header(self,title=None,menu=None):
        """show headings and column names for DB"""

        imgdir = self.imgdir
        html_header = 'Content-Type: text/html; charset=utf-8\n'
        html_header += '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/loose.dtd">\n\n'
        print html_header
        print '<html>'
        print '<head>'

        if title==None:
            print '<title>PEAT Web Interface</title>'
        else:
            print '<title>%s</title>' % title
        print '<link href="%s/styles.css" rel="stylesheet" type="text/css" />' %self.bindir
        print '<link rel="shortcut icon" href="%s/favicon.ico" type="image/x-icon" />' %self.bindir
        print '<script type="text/javascript" src="/scripts/checkbox.js"></script>'
        print '</head>'
        print '<body>'
        print '<div class="header">' 
        print '<p id="title">PEAT Web Interface: %s</p>' %self.project
        print '</div>'
        print '<script type="text/javascript" src="/scripts/boxover.js"></script>'
        print '<hr>'
        print '<div>'

        if menu==1:
            self.menu()

        return

    def write_sessionkey(self, action=None, project=None, manage=None):
        """Write the session key"""
        if self.user:
            print '<input type=hidden name=login value=%s>' %self.user
            print '<input type=hidden name=sessionkey value=%s>' %self.sessionkey
        if action:
            print '<input type=hidden name=action value=%s>' %action
        if project:
            print '<input type=hidden name=project value=%s>' %project
        if manage:
            print '<input type=hidden name=manage value=yes>'

        return

    def footer(self):
        """Print page footer"""

        print '<br>'
        print '<p><center>Provided by <a href="http://enzyme.ucd.ie">Nielsen Research Group at UCD</a></center></p>'
        print '</div>'
        print '</div>'
        print '</body>'
        return

    def menu(self):
        """Print the menu"""
        bindir = self.bindir
        print '<div class="menu">'
        print '<table id="menu" valign=top align=left>'
        print '<th><b>Menu</b></th>'

        print '<form action="%s/main.cgi" METHOD="POST" ENCTYPE="multipart/form-data">' %bindir
        self.write_sessionkey('show_intro')
        print '<tr><td class="menu">\
            <input type=submit value="Home" name=submit size=20 class="btn"></td></tr>'
        print '</form>'
        print
        print '<form action="%s/main.cgi" METHOD="POST" ENCTYPE="multipart/form-data">'  %bindir
        self.write_sessionkey('show_all')
        print '<tr><td class="menu">'
        print '<input type=submit value="Show All Records" name=submit class="btn"'
        print """title="header=[Show all] body=[Show all records in the DB]"></td></tr>"""
        print '</form>'

        print
        print '<form action="%s/main.cgi" METHOD="POST" ENCTYPE="multipart/form-data">' %bindir
        self.write_sessionkey('summary')
        print '<tr><td class="menu">\
            <input type=submit value="Summary" name=submit class="btn"></td></tr>'
        print '</form>'

        print '<form action="%s/main.cgi" METHOD="POST" ENCTYPE="multipart/form-data">' %bindir
        self.write_sessionkey('show_help')
        print '<tr><td class="menu">\
            <input type=submit value="Help" name=submit class="btn"></td></tr>'
        print '</form>'
        print

        searchmessage = 'Enter your search here'
        print '<tr><td class="menu"><a><b>SEARCH</a></td></tr>'
        print '<form action="%s/main.cgi" METHOD="POST" ENCTYPE="multipart/form-data">' %bindir
        self.write_sessionkey('search')
        #protein search box
        print '<tr><td class="menu">'
        print 'Protein: '
        print '<select name="proteinmatchmethod" class="btn">'
        print '<option value="or">Any words'
        print '<option value="and">All words'
        print '</select>'
        print '</td></tr>'
        print '<td class="menu">'
        print """<input type="text" name="words" size=22 maxlength=50 value="" \
                title="header=[Protein search] body=[Enter multiple names seperated by a space.\
                    This search is not case sensitive]"></td></tr>"""

        #dataset search box
        print '<tr><td class="menu">'
        print 'Dataset: '
        print '<tr><td class="menu">'
        print """<input type="text" name="dataset" size=22 maxlength=40 value="" class="btn"\
                 title="header=[Name of data item] body=[Use any strings]"></td></tr>"""


        #drop list for whether to match any or all of the searches from each box
        print '<tr><td class="menu">'
        print """Match: <select name="matchmethod" class="btn" title="header=[Global match method] \
                 body=[Match records with ALL of the search criteria or only those with ANY those\
                  attributes]">"""
        print '<option value="and">All'
        print '<option value="or">Any'
        print '</select> <p>'
        print '</td></tr>'
        print '</td></tr>'
        print '<tr><td class="menu"><input type="submit" value="Search Now" class="btn"></td></tr>'
        print '</form>'
        print '</table>'
        print '</div>'
        print '<div id="content">'
        return

    def show_all(self):
        '''Show all records in DB
           From the locally saved project if present, otherwise checkout'''

        self.DB = self.connect()
        self.show_DB_header(menu=1)
        sys.stdout.flush()
        self.show_DB()
        return

    def connect(self):
        """Connect to the peat db"""        
        saveout = sys.stdout
        fsock = open('out.log', 'w')
        sys.stdout = fsock
        sys.stdout.flush()
        print self.server
        DB = PDatabase(server=self.server, port=self.port,
                      username=self.user,
                      password=self.password,
                      project=self.project)   
        sys.stdout.flush()
        sys.stdout = saveout
        fsock.close()     
        return DB
        

    def show_DB(self, selected=None):
        """
        Show all the records in the PEAT database in a html table.
        """
        DB = self.DB
        ptmodel = PEATTableModel(DB)
        colnames = ptmodel.columnNames
        exclude = ['pdb Info','Mutations']
        ekincols = ptmodel.ekintypes
        print '<table id="mytable" cellspacing=0>'
        print '<tr>'
        #do column names at top first
        for c in colnames:
            if c in exclude:
                continue
            print '<th>'+c+'</th>'
        print '</tr>'

        y_count=1
        r=0
        #used for summarising titration curve
        tit1H=0
        tit15N=0
        tit13C=0
        #get sort order by name field
        if selected == None:
            selected = DB.getRecs()
        names = []
        for p in selected:
            names.append((DB[p].name, p))
        names.sort()

        for name in names:
            protein = name[1]
            if r % 2 == 0:
                cls = 'alt'
            else:
                cls = ''
            print '<tr>'

            for column in colnames:
                if column in exclude:
                    continue
                #get field type here
                fieldtype=None
                if DB['userfields'].has_key(column):
                    fieldtype=DB['userfields'][column]['field_type']
                   
                elif column == 'Structure':
                    fieldtype='structure'
                if fieldtype==None:
                    fieldtype = 'text'

                display_text = DB[protein].getDisplayAttribute(column)
                if display_text == '':
                    display_text = '-'
                if fieldtype == 'Text' or fieldtype == 'text' or fieldtype == 'structure':

                    if column == 'Name':
                        print '<th class="spec"> %s </th>' % display_text
                    else:
                        print '<td class=%s> %s </td>' % (cls,display_text)
                #links
                elif fieldtype == 'Link':
                    display_text = '-'
                    if DB[protein].has_key(column) and  isinstance(DB[protein][column], dict):
                        if DB[protein][column].has_key('text'):
                            display_text = DB[protein][column]['text']
                            link_text = DB[protein][column]['link']
                            if link_text != None and link_text != '':
                                if column == 'PDB_link':
                                    viewlink = 'http://firstglance.jmol.org/fg.htm?mol='+display_text
                                    print '<td class=%s> <a href=%s target="_blank">%s</a> <a href=%s target="_blank"><img class="thumb" \
                                        src="%s/fg.png"></a></td>' % (cls,link_text,display_text,viewlink, self.imgdir) 
                                elif column == 'PMID_link':
                                    authors=''
                                    try:
                                        auth, title = DB[protein][column]['authors'], DB[protein][column]['title']
                                        for a in auth:
                                            authors+=a.encode("utf-8")+' '
                                    except:
                                        auth, title = ('NA', 'NA')
                                    print '<td class=%s> <a href=%s target="_blank" title="header=[%s] \
                                            body=[%s]">%s</a></td>' % (cls,link_text,authors,title,display_text)
                                else:
                                    print '<td class=%s> <a href=%s target="_blank">%s</a></td>' % (cls,link_text,display_text)
                            else:
                                print '<td class=%s> %s </td>' % (cls,display_text)
                    else:
                        print '<td class=%s> %s </td>' % (cls,display_text)
                #ekin data fields
                elif fieldtype in ekincols:
                    if display_text.startswith('  0 recs') or fieldtype == 'Text':
                        print '<td class=%s> %s </td>' % (cls, display_text)
                    else:
                        import urllib
                        urlfields = urllib.urlencode({'login':self.user,'project':self.project,
                                           'sessionkey':self.sessionkey,
                                           'protein':protein,'column':column,
                                           'fieldtype':fieldtype,'action':'show_specific'})
                        print '<td class=%s> <a href="%s/main.cgi?%s" target="_blank">%s</a></td>' % (cls,self.bindir,urlfields,display_text)
                        #get no. of curves in that record
                        try:
                            numcurves = int(display_text.split('recs')[0])
                        except:
                            numcurves = 1
                        if '1H' in column:
                            tit1H=tit1H+1*numcurves
                        elif '15N' in column:
                            tit15N=tit15N+1*numcurves
                        elif '13C' in column:
                            tit13C=tit13C+1*numcurves

            print '</tr>'
            r=r+1
            # Update the row count
            y_count=y_count+1

        print '<tr><td valign=top colspan=3> summary: %s proteins 1H=%s 15N=%s 13C=%s</tr>' %(r, tit1H, tit15N, tit13C)
        print '</table>'
        self.footer()
        return


    def show_specific(self):
        """Show datasets for an individual cell, depends on fieldtype"""

        protein = self.form.getfirst('protein')
        rectype = self.form.getfirst('rectype')
        col = self.form.getfirst('column')
        self.show_DB_header(title=str(protein+': '+col), menu=1)
        sys.stdout.flush()

        DB = self.DB = self.connect()    
        ptmodel = PEATTableModel(DB)
        ekincols = ptmodel.ekintypes
        fieldtype = DB['userfields'][col]['field_type']
        protname = DB[protein]['name']
        print '<div><a>Protein: %s </a><br>' %protname
        print 'Column data: %s </div>' %col
        from PEATDB.Ekin.Web import EkinWeb

        if fieldtype in ekincols:
            fieldtype = 'ekindata'
        #plot ekin data curves, if fieldtype is ekindata
        if fieldtype == 'ekindata':
            E = DB[protein][col]           
            EW = EkinWeb()
            EW.showEkinPlots(project=E, datasets='ALL', path=self.plotsdir,                            
                             imgpath=self.imagepath)
            print '</table>'
        #save the pdb structure to a file and provide url link to it
        elif fieldtype == 'structure':
            recdata = DB[protein][col]
            for k in recdata:
                print '<br><a> %s</a>' %str(k)
        elif fieldtype == 'Notes':
            recdata = DB[protein][col]
            print '<br><a> %s</a>' %recdata['text']      
        elif fieldtype == 'File':
            recdata = DB[protein][col]
            print '<p> This record contains a data file. Click link to open.'
            print '<a href="http://peat.ucd.ie/titration_db/%s">Open link</a>' %recdata['location']
        else:
            print '<br>Sorry, no handler for this data right now..</br>'
        self.footer()
        return

    def show_datasets(self):
        """Show single or multiple dataset plot from search results"""
        plotoption = self.form.getfirst('plotoption')
        norm = bool(self.form.getfirst('normalise'))
        logx = bool(self.form.getfirst('logx'))
        logy = bool(self.form.getfirst('logy'))
        legend = bool(self.form.getfirst('legend'))
       
        DB = self.DB = self.connect()       
        from PEATDB.Ekin.Web import EkinWeb
        ekincols = DB.ekintypes

        self.show_DB_header(menu=1)
        self.menu()
        sys.stdout.flush()

        #get dataset name(s) from form
        residuefield = self.form.getvalue('residue')

        if not isinstance(residuefield, list):
            residuefield = [residuefield]

        #prepare data, if more than one dataset from seperate proteins,
        print '</div>'
        print '<table id="mytable" cellspacing=0>'
        print '<tr><th>'
        print 'plotoption:', plotoption
        print 'norm:', str(norm)
        print 'log-x', logx
        print 'log-y', logy
        print 'legend', legend
        print '</table>'

        ekindata={}
        ekindata['__datatabs_fits__']={}
        ekindata['__meta_data__']={}
        #extracting seperate datasets and put them in ekindata
        for r in residuefield:
            fields = r.split('$')
            protein = fields[0]
            col = fields[1]
            d = fields[2]
            E = DB[protein][col]
            fits = E.__datatabs_fits__       
            meta = E.__meta_data__
            protname = DB[protein]['name']
            ekindata['__datatabs_fits__'][d+'_'+col+'_'+protname] = fits[d]
            ekindata['__meta_data__'][d+'_'+col+'_'+protname] = meta[d]
            ekindata[d+'_'+col+'_'+protname] = E.getDataset(d)

        if plotoption == None:
            plotoption = 1
        else:
            plotoption = 3

        EW = EkinWeb()
        EW.showEkinPlots(ekindata, 'ALL', path=self.plotsdir,
                        imgpath=self.imagepath,
                        plotoption=plotoption, normalise=norm, legend=legend,
                        logx=logx, logy=logy)

        return

    def do_search(self, globalop, proteinop, proteins, datasets):
        """Do searches for the various parameters, name, pka etc"""

        import re, urllib
        from PEATDB.PEATTables import PEATTableModel
        from PEATDB.Ekin.Fitting import Fitting
        from PEATDB.Ekin.Web import EkinWeb
        DB = self.DB  
        EW = EkinWeb()
        ekincols = DB.ekintypes

        found=[]
        protlist = proteins.split(' ')
        residuelist = residues.split(' ')

        def createPhrase(items, op):
            if op == 'or':
                logicsymbol = '|'
            else:
                logicsymbol = '+'
            itemlist = items.split(' ')
            phrase =''
            c=1
            for i in itemlist:
                if c == len(itemlist):
                    phrase = phrase + i
                else:
                    phrase = phrase + i + logicsymbol
                c=c+1
            return re.compile(phrase, re.IGNORECASE)

        #if there is a protein name entered, get the list of proteins/records to use
        #otherwise use all records to search for the other keys
        names = []
        keywords = {}
        for p in DB.getRecs():
            name = DB[p].name
            names.append((name,p))
            if hasattr(DB[p], 'keywords'):
                keywords[name] = DB[p].keywords
            else:
                keywords[name] = ''
        names.sort()

        #create search expressions
        s_protein = createPhrase(proteins, proteinop)
        s_residue = createPhrase(residues, 'or')
        pkarange = pka.split('-')

        #do search
        for name in names:
            proteinname = name[0]
            protein = name[1]
            if s_protein.search(proteinname) or s_protein.search(keywords[proteinname]):
                found.append(protein)

        if len(found)==0:
            print '<h2>Sorry, no proteins found that match this name. <h2>'
            return
        elif residues == '' and pka == '':
            self.show_DB(selected=found)
            return found

        #now search for the other keys inside selected proteins if needed
        ekinfound={}; foundfits={}
        ignore =['__fit_matches__', '__Exp_Meta_Dat__', '__datatabs_fits__']
        
        kys = list(DB.getRecs())
        kys.sort()
        if len(found) == 0 or globalop == 'or':
            found = DB.getRecs()

        print '<table id="mytable" cellspacing=0 align=center>'
        '''search recs for residue and pkarange
           and add dataset names to new ekinfound dict'''
        for protein in found:
            for col in DB[protein].keys():
                if nucleus != 'any':
                    if nucleus not in col:
                        continue
                E = DB[protein][col]
                if DB['userfields'].has_key(col):
                    fieldtype=DB['userfields'][col]['field_type']
                else:
                    fieldtype=None
                if fieldtype in ekincols:
                    dk = E.datasets                    
                    fits = E.__datatabs_fits__                    
                    for k in dk:
                        meta = E.getMetaData(k)
                        try:
                            thisres = meta['residue']
                        except:
                            thisres = k         
                        if thisres == None:
                            thisres = k
                        if residues != '':
                            if s_residue.search(thisres) and not k in ignore:
                                if not ekinfound.has_key(protein+'$'+col):
                                    ekinfound[protein+'$'+col]=[]
                                ekinfound[protein+'$'+col].append(k)

                        if pka != '':
                            foundpka = 0
                            if k in fits.keys():
                                if fits[k] == None:
                                    continue
                                if fits[k].has_key('model'):
                                    model = fits[k]['model']
                                else:
                                    continue
                                if not foundfits.has_key(protein+'$'+col):
                                    foundfits[protein+'$'+col]={}
                                X = Fitting.getFitter(model)
                                pnames = X.names
                                i=0
                                #check if first num is larger, then swap
                                try:
                                    pk1=float(pkarange[0])
                                    pk2=float(pkarange[1])
                                except:
                                    print '<h2> Error: pka values are not valid </h2>'
                                    return
                                if pk1 > pk2:
                                    tmp = pk1
                                    pk1 = pk2
                                    pk2 = tmp

                                #iterate thru parameter names and match any pK fields
                                #this code is not that efficient!
                                for p in pnames:
                                    if 'pK' in p:
                                        pka = fits[k][i]
                                        if pka >= pk1 and pka <= pk2:
                                            foundpka=1
                                    i=i+1
                                #if match is 'ANY', just append dataset if not there already
                                #also for case if no residue value entered
                                if globalop == 'or' or residues == '':
                                    if foundpka == 1:
                                        if not ekinfound.has_key(protein+'$'+col):
                                            ekinfound[protein+'$'+col]=[]
                                        if not k in ekinfound[protein+'$'+col]:
                                            ekinfound[protein+'$'+col].append(k)
                                #if match is 'ALL', need to check dataset already found
                                #and if no pka found, remove it. if both are true keep it
                                elif globalop == 'and':
                                    if foundpka == 0:
                                        if ekinfound.has_key(protein+'$'+col) and k in ekinfound[protein+'$'+col]:
                                            #print 'removing', protein, col
                                            ekinfound[protein+'$'+col].remove(k)

                                foundfits[protein+'$'+col][k]=fits[k]
        #check for empty fields in ekinfound dict and delete them..
        for d in ekinfound.keys():
            if len(ekinfound[d])==0:
                del(ekinfound[d])

        #if no results, just say that and return
        if len(ekinfound)==0:
            print '<h2> Sorry, no records found that match these parameters. </h2>'
            return
        #top display button and options
        print '<form name="resultsform" action="%s/main.cgi" METHOD="POST" ENCTYPE="multipart/form-data">' %self.bindir
        self.write_sessionkey('show_datasets')
        print '<td valign=top align=right colspan=3> <input type=submit value="display selected" name=submit>'
        print '<label><input type=checkbox name="plotoption" value=3>single graph</label>'
        print '<label><input type=checkbox name="normalise" value=1>normalise</label>'
        print '<label><input type=checkbox name="logx" value=1>log-x</label>'
        print '<label><input type=checkbox name="logy" value=1>log-y</label>'
        print '<label><input type=checkbox name="legend" value=1>legend</label>'
        print '<input type=button value="Check All" onClick="checkAll(document.resultsform.residue)">'
        print '<input type=button value="Uncheck All" onClick="uncheckAll(document.resultsform.residue)"></td>'
        print '</tr>'
        #header
        cols=['protein','column','residues']
        for c in cols:
            print '<th>'+c+'</th>'
        print '</tr>'
        ekys = ekinfound.keys()
        ekys.sort()

        r=1
        for k in ekys:
            if r % 2 == 0:
                cls = "spec"
            else:
                cls = ""

            fields = k.split('$')
            protein = fields[0]
            column = fields[1]

            proteinname = DB[protein]['name']
            print '<tr>'
            print '<th class="spec">%s</th> <td> %s </td>' %(proteinname, column)
            print '<td>'
            for residue in ekinfound[k]:
                residuefield = k+"$"+residue
                fithtml = ''
                try:
                    print '<input type=checkbox id="residue" name="residue" value=%s>' %("'"+residuefield+"'")
                except:
                    print 'UnicodeEncodeError'
                    continue

                urlfields = urllib.urlencode({'login':self.user,'project':self.project,
                                           'residue': residuefield,'action': 'show_datasets'})
                print '<a href="/cgi-bin/titration_db/main.cgi?%s" target="_blank">%s</a>'\
                          % (urlfields, residue)

            print '</tr>'
            r=r+1
        print '</form>'
        print '</table>'
        self.footer()

        return

    def show_search_results(self):
        """Display search results"""

        self.show_DB_header(menu=1)
        sys.stdout.flush()
        self.DB = self.connect()
       
        key=self.form.getfirst('key')
        matchmethod=self.form.getfirst('matchmethod')
        proteinmatchmethod=self.form.getfirst('proteinmatchmethod')
        words=self.form.getfirst('words')
        residue=self.form.getfirst('residue')
        nucleus=self.form.getfirst('nucleus')
        pka=self.form.getfirst('pka')

        print '<table bgcolor=#CD9B9B border="1" bordercolor=black cellspacing=0 cellpadding=4 valign=top \
                align=center width=80%><tr><td>'

        print '<div align=left>'
        print '<big><a>Search results for protein=%s %s residue=%s %s nucleus=%s pka=%s</big></a>'\
                  %(words,matchmethod,residue,matchmethod,nucleus,pka)
        print '<br>'
        print '</div>'
        print '</table>'
        self.do_search(matchmethod, proteinmatchmethod, words, residue, nucleus, pka)
        return


if __name__ == '__main__':
    #test
    PEATWeb(project='titration_db',
          user='guest',
          passwd='123',
          bindir='/cgi-bin/newtitdb',
          fullpath='/var/www/html/newtitdb')
    

