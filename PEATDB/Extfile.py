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


"""
Classes to handle external files via peatdb
Author: Damien Farrell 2010
"""

from Tkinter import *
from GUI_helper import *
import os, types
import tkSimpleDialog, tkFileDialog, tkMessageBox
import platform
import Pmw
from Record import FileRecord


class FileHandler(object):
    """Class to handle external files for peatdb, uses some stuff from
       old peat from db_actions etc
       Should add thumbnail stuff here??"""

    extensions = {'image' : ['.jpg','.jpeg','.tiff','.tif','.gif','.ps','.png','.ico'],
                  'doc' : ['.doc','.xls','.ppt','.pdf','.docx','.txt','.csv','.odt'],
                  'compressed' : ['.zip','.gz','.bz','.tar'],
                  'any':''}

    def __init__(self, parent=None, DB=None):
        self.parent = parent
        self.currplatform = platform.system()

        return

    def displayFile(self, filerec):
        """Uses a new method of file storage with blobs"""
        #now we just get the blob from the filerecord object...
        from PILView import PILViewer
        from PEATDB.Record import FileRecord
        #check class        
        if not type(filerec) is FileRecord:           
            return False
        myblob = filerec.blob
        ext = filerec.ext
        mtype = filerec.mimetype
        print 'mimetype is', mtype
        if myblob != None:
            f = myblob.open("r")
            filename = f.name
            f.close()
        else:
            return False

        if self.currplatform == 'Linux':
            #try using mailcaps on linux system
            import mailcap, os
            d=mailcap.getcaps()
            print mtype, filename
            f=mailcap.findmatch(d, mtype)[1]
            os.system(f['view'] %filename)
            return
        else:
            import os, tempfile, shutil
            tempname = os.path.normpath(os.path.join(tempfile.gettempdir(), filerec.name))
            if not os.path.exists(tempname):
                #os.remove(tempname)
                shutil.copy(filename, tempname)
            os.startfile(tempname)

        return True

    def addFileDialog(self, filerec=None):
        """Do the front-end stuff for adding files to peat"""
        import tkFileDialog
        filename=tkFileDialog.askopenfilename(defaultextension='.*',
                                              initialdir=os.getcwd(),
                                              filetypes=[("All files","*.*")])
        self.newfile = filename
        return filename


    def importFileset(self, DB, parentframe=None, callback=None):
        """Import a series of external files in a folder to peat"""

        def getdir():
            import tkFileDialog, os
            d=tkFileDialog.askdirectory(initialdir=importdir.get(),title='Select directory files',
                                                mustexist=1,parent=self.parent.master)
            if d:
                importdir.set(d)
                return

        def doimport():
            idir = importdir.get()
            ext = filetype.get()
            if idir == None:
                return
            files = os.listdir(idir)
            okfiles = []
            for f in files:
                pth = os.path.join(idir,f)
                if os.path.isfile(pth):
                    if ext != 'any':
                        if os.path.splitext(pth)[1] not in self.extensions[ext]:
                            continue
                    okfiles.append(pth)
            print 'files to be used:', okfiles
            DB.addField('file', 'File')
            from Dialogs import PEATDialog
            self.pb=PEATDialog(self.parent.master, option='progressbar',
                                          message='Importing files')
            self.pb.update_progress(0)
            total = len(okfiles)
            row=0
            #iterate over filenames in folder and add one rec to DB for each file
            for f in okfiles:
                recname=os.path.basename(f)
                DB.addBlob(recname, 'file', f)
                print 'added %s as %s' %(f,recname)
                row=row+1
                c=float(row)/float(total)*100.0
                self.pb.update_progress(c)
            self.pb.close()
            if callback!=None:
                callback()
            return

        def close():
            win.destroy()
            return

        if parentframe!=None:
            win = Frame(master=parentframe,relief=RAISED)
            win.pack(fill=BOTH)
        else:
            win = Toplevel()
            win.title('Import ext files')

        importdir=StringVar(); importdir.set(os.getcwd())
        Label(win,textvariable=importdir,bg='white',justify=RIGHT).pack(fill=BOTH,side=TOP)
        Button(win,text='Select Directory',command=getdir).pack(fill=BOTH,side=TOP)
        filetype=StringVar(); filetype.set('any')
        optmenu = Pmw.OptionMenu (win,
                labelpos = 'w',
                label_text = 'File Type:',
                menubutton_textvariable = filetype,
                items = self.extensions.keys(),
                menubutton_width = 10,
                #command = ,
        )
        optmenu.pack(fill=BOTH,side=TOP)
        Button(win,text='Import',command=doimport).pack(fill=BOTH,side=TOP)
        Button(win,text='Cancel',command=close).pack(fill=BOTH,side=TOP)

        return

    def exportExtFiles(self, DB):
        """Get all stored ext files in DB and save to folder"""
        import tkFileDialog, os, shutil
        expdir=tkFileDialog.askdirectory(initialdir=os.getcwd(),title='Select a directory',
                                          mustexist=1,parent=self.parent.master)
        if not expdir:
            return
        path = os.path.join(expdir,'peat-export')
        if not os.path.exists(path):
            os.mkdir(path)
        count=0
        files = DB.getExtFilenames()
        for f in files:
            shutil.copyfile(files[f], os.path.join(path, f))
            count+=1
        tkMessageBox.showinfo("Exported",
                             'Done. Exported %s files to\n%s'%(count,path))
        return

    def displayExtFiles(self, DB):
        """Display external files (blobs) in DB"""
        files = DB.getExtFilenames()
        print files
        if len(files) == 0:
            return
        items=[]
        for f in files:
            print f, files[f]
            items.append('<a href=%s>%s</a>' %(files[f],f)
                         )
        import markup
        page = markup.page( )
        page.init( title="PEAT_DB Files Preview",
                   css=( 'one.css' ),
                   header="Preview for project",
                   footer="PEAT_DB" )
        page.ul( class_='mylist' )
        page.li( items, class_='myitem' )
        page.ul.close()

        filename = '/tmp/prevtemp.html'
        hf = open(filename,'w')
        for c in page.content:
            hf.write(c)
        hf.close()
        import webbrowser
        webbrowser.open(filename, autoraise=1)
        return


