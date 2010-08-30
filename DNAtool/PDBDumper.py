#!/bin/env python
#
# Class for dumping text files of current primer DB in DNAtool
#
# (C) Copyright 2005-2006 Jens Erik Nielsen, University College Dublin. All rights reserved
#

import sys, os,copy
from Tkinter import *

class PDBDumper:

    def __init__(self,parent,Dump2Dir,noDB=None):
        from datetime import datetime
        self.MaxDumpFiles = 3
        self.parent = parent
        #project name from DB_Main.. too many levels of reference!!
        #should set the proj name from where this objected is created..
        if noDB==None:
            projName = self.parent.parent.parent.DB.meta['info']['project']
        else:
            projName = 'current'
        now = datetime.today()
        print " projname ",projName,"\n\n",str(now.time())," ",str(now.time()).split('.'),now.date(),
        self.DumpName = projName+'.'+"PrimerDB"+"."+str(now.date())+str(now.time()).split('.')[0]+".csv"
        self.DumpBase = projName+'.'+"PrimerDB"
        print " DumpName "+self.DumpName
        self.DIR = Dump2Dir
        self.CheckDump()
        #current primer db held in the parent DNAtool instance
        if noDB==None:
            self.primerDB = self.parent.parent.data['primer_dict']
        else:
            self.primerDB = []
        return

    # Does the file writing
    def doDump(self,event=None):
        """Does the file writing"""
        # if pdb is empty do nothing
        try:
            self.parent.parent.data.get('primer_dict')
        except:
            return
        if event.widget == self.parent.pDB_win:
            print "Doing dump.."
            self.CheckDump()

            list_csv = self.createCSV(self.primerDB)
            #DFILE = self.DIR+'/'+self.DumpName
            DFILE =os.path.join(self.DIR,self.DumpName)
            #print "writing: ",list_csv
            try:
                fd=open(DFILE,'w')
                for line in list_csv:
                    str = ''
                    for S in line:
                        str = str+S+","
                    fd.write(str+'\n')
                fd.close()
                print "File written ",DFILE
            except:
                print "could not write dump or pdb was empty"
        # open file

        #write file

        return

    # Creates the CSV file structure
    def createCSV(self,Labbook):
        """Creates the CSV file structure"""
        result = None
        if self.parent.parent.data.has_key('primer_dict'):
            PDB = self.parent.parent.data['primer_dict']
        else:
            return
        primer_names=PDB.keys()
        HeaderList = ['description', 'sequence']
        ActualHeaderList = ['name','description', 'sequence']
        DumpPDB = []
        #DumpPDB.append(ActualHeaderList)
        #print "Header list ",HeaderList
        for entry in PDB:
            #print "entry ",entry
            tmp = []
            tmp.append(entry)
            for head in HeaderList:
                #print PDB[entry][head]
                tmp.append(PDB[entry][head])
            DumpPDB.append(tmp)

        result = DumpPDB
        #print " DumpPDB ",DumpPDB
        return result

    # Checks to see which is oldest version present and overwrites it
    def CheckDump(self):
        import os,os.path
        print "Checking Dump"
        Files = os.listdir(self.DIR)
        N_files = 0
        Dump_files = []
        for count,I in enumerate(Files):
            if I.find(self.DumpBase) != -1:
                #print "file ",self.DIR+"/"+I
                Dump_files.append([I,os.path.getmtime(self.DIR+"/"+I)])
                N_files += 1
        Dump_files.sort(lambda x, y: cmp(x[1], y[1]))
        Dump_files.reverse()
        #print "Dumps ",Dump_files
        # if we have more dump files than MaxDumpFiles, delete the oldest
        # and create a new one
        try:
            #print "N ",len(Dump_files)
            while len(Dump_files) >= self.MaxDumpFiles:
                del_file = Dump_files.pop()
                #print "remove dump file :",self.DIR+"/"+del_file[0],"\n\n"
                rfile=os.path.join(self.DIR,del_file[0])
                os.remove(rfile)
        except:
            print "Could not remove dump file "
           #Dump_files.sort( os.getmtime(  	path) )

        #getmtime(  	path)
        return

