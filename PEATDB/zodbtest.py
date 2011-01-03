#!/usr/bin/env python
#
# This file is part Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-
# All rights reserved
#
# Author: Damien Farrell 2009

#testing for ZODB/ZEO as backend PEAT replacement

import os, sys, random
import numpy
import timeit
import ZODB, transaction
from ZODB import FileStorage, DB
from BTrees.OOBTree import OOBTree
from BTrees.OOBTree import OOBucket
from ZEO.ClientStorage import ClientStorage
from ZODB.blob import Blob
from persistent import Persistent
from ZODB.PersistentMapping import PersistentMapping
from ZODB.PersistentList import PersistentList
from PEATDB.Base import PDatabase, zDatabase
from PEATDB.Record import PEATRecord, Meta, FileRecord
from PEATDB.Ekin.Base import EkinProject
from Actions import DBActions
from time import time

ekintypes = ['General', 'NMR titration', 'Protein Stability',
                      'pH Activity profile', 'Michaelis-Menten kinetics',
                      'pH stability profile', 'Amyloid formation']


class Item(Persistent):
    """Item inside a record has to be put into this class so we can
       have lazy loading of keys in a record"""
    def __init__(self, data):
        self.data = data
        return
        
class NewRecord(OOBTree):
    """Proto class to represent a persistent object in PEAT.
       We subclass this for protein records and meta info"""
    def __init__(self, name):
        OOBTree.__init__(self)
        self.name=name
        return

    def getDisplayAttribute(self, key, fieldtype):
        """Get the text to display for this attribute of the record
           This is used for the table display to cache the values"""
        text = ''
        if not hasattr(self, key):
            return ''
        data = self[key]
        if data == None:
            return ''
        if fieldtype == 'text':
            text = data
        return text

    def getFields(self):
        return self.keys()

    def addField(self, key, data=None):
        self[key] = data
        return

    def delete(self, key):
        del self[key]
        self._p_changed = 1
        return

    def getField(self, key):
        return self[key]

    def get_type(self, key):
        return type(self[key])
        raise KeyError('No such key: %s' %key)

    def __getitem__(self, key): 
        return OOBTree.__getitem__(self, key).data

    def __setitem__(self, key, data):        
        item = Item(data)
        OOBTree.__setitem__(self, key, item)       
        return
        
    def __getattr__(self, key):
        """Maps values to attributes.
        Only called if there *isn't* an attribute with this name
        """
        if Persistent._p_getattr(self, key):
            return Persistent.__getattribute__(self, key) 
        try:
            return self.__getitem__(key)
        except KeyError:
            raise AttributeError(key)

    def __setattr__(self, key, data):
        """Maps attributes to key values."""
        if Persistent._p_setattr(self, key, data):
            return
        self._p_changed = 1
        return self.__setitem__(key, data)
        

def simulate_sequence(length):
    """simulate a dna seq"""
    dna = ['A', 'C', 'G', 'T']
    sequence = ''
    for i in range(length):
        sequence += random.choice(dna)
    return sequence

def Timeit(func,number=1,module="__main__"):
    """ A wrapper which can be used to time any function """

    name = func.__name__
    t = timeit.Timer("%s()"%name, "from %s import %s" % (module, name))
    return "%.2f usec/pass" % (1000000*t.timeit(number=number)/number)

def speedTest1():
    """benchmark read and write"""
    DB = PDatabase(server='localhost',port=8090)
    #DB = PDatabase(local='./Data.fs')
    for i in range(100,300):
        DB.add(i)
    for i in range(100,290):
        DB.delete(i)
    DB.commit()
    print DB
    DB.close()
    return

def createdb(local=None, server=None, project=None, username=None, norecs=1000):
    """Create and add some test data"""
    if local != None:
        if os.path.exists(local):
            for i in ['.lock','.index','']:
                try:
                    os.remove(local+i)
                except:
                    pass
        DB = PDatabase(local=local)
    elif server!=None:
        DB = PDatabase(server=server, username=username,
                    password='123', project=project)

    import string
    import DNAtool.mutation as mutation
    
    choices = ['a','b','c','d']
    DB.addField('choice', 'text')
    DB.addField('stab', 'text')
    DB.addField('activity', 'text')
    #DB.addField('ekin', 'General')
    E = EkinProject()
    data=E.readDataset('Ekin/test.Ekindat')
    E.insertDataset(data['data'], 'test', fit=data['fit'])
    print 'creating dummy data..'
    j=0
    count=0
    for i in range(norecs):
        if j>3: j=0
        c=''
        for k in range(10):
            c += random.choice(string.letters)
        DB.add(c)
        DB.data[c].choice = choices[j]
        DB.data[c].DNASeq = simulate_sequence(300)
        AAseqs3,AAseqs1 = mutation.translate(DB.data[c].DNASeq)
        DB.addProtseq(c, AAseqs3[1][5:80], 1)
        DB.data[c].stab = str(round(random.normalvariate(1,2),3))
        DB.data[c].activity = str(round(random.normalvariate(30,4),3))
        #DB.data[c].ekin = E
        j+=1
        count+=1
        if count>3000:
            print 'saving..'
            DB.commit()
            DB.db.cacheMinimize()
            count=0

    DB.commit()
    return DB

def checkdb(server=None, project=None, local=None):
    """test basic ops for the db"""
    if server != None:
        DB = PDatabase(server=server, port=8080,
                  username='farrell',
                  password='123',
                  project=project)
    else:
        DB = PDatabase(local=local)

    print DB
    #DB.listRecs()
    '''first = DB.data.keys()[0]
    rec = DB.data[first]    
    for k in rec.keys():
        print type(rec[k])'''
    print DB['userfields']
    print DB.meta

    DB.close()
    return

def importOldProj(datadir,local=None, server=None,
                    project=None, username=None):
    """Import old peat projects"""
    import PEAT_DB.Database as peatDB
    from PEAT_DB.PEAT_dict import PEAT_dict, sub_dict    
    import copy
    if local != None:
        newDB = PDatabase(local=local)
    elif server != None:
        newDB = PDatabase(server=server, username=username, port=8080,
                          password='123', project=project)

    print newDB
    PT = peatDB.Database(datadir, Tk=False)
    oldDB = PT.DB
    print 'got old peat_db with %s proteins' %len(PT.proteins)

    print PT.DB.keys()
    #import meta stuff like userfields, table
    for p in newDB.meta.special:
        if not p in PT.DB.keys():
            continue
        print 'adding',p
        for k in PT.DB[p]:
            newDB.meta[p][k] = copy.deepcopy(PT.DB[p][k])
    newDB.meta._p_changed = 1

    for p in PT.proteins:
        if p in newDB.meta.special:
            continue

        name = oldDB[p]['Name']         
        rec = PEATRecord(name=name)
        for col in oldDB[p].keys():
            cdata = oldDB[p][col]
            recdata = {}
            if col == 'name':
                cdata = oldDB[p]['Name']
  
            if oldDB['userfields'].has_key(col) and oldDB['userfields'][col]['field_type'] in ekintypes:
                E=EkinProject(data=cdata)
                E.length = len(E.datasets)
                if len(E.datasets)==0:
                    continue
                cdata = E

            if type(cdata) == sub_dict:
                for k in cdata.keys():
                    recdata[k] = copy.deepcopy(cdata[k])
            else:
                recdata = cdata
            if cdata != '' and cdata != None:
                rec.addField(col, data=recdata)
        newDB.add(p,rec)
    print newDB.meta.userfields
    #remove any file cols, too hard to import
    for m in newDB.meta.userfields.keys()[:]:
        if newDB.meta.userfields[m]['field_type'] == 'File':
            newDB.deleteField(m)
    newDB.commit(user='farrell', note='import')
    newDB.close()
    print 'import done'

    return

def setDisplayFields(local=None, server=None,
                    project=None, username=None):
    """Update display attributes for every cell"""
  
    if local != None:
        DB = PDatabase(local=local)
    elif server != None:
        DB = PDatabase(server=server, username=username, port=8080,
                          password='123', project=project)
        
    for r in DB.getRecs():
        rec = DB[r]
        for f in rec:
            x = rec.setDisplayAttribute(f, rec[f])
            rec._display[f] = x
            print f, rec.getDisplayAttribute(f)
    #DB.commit('displ attrs')        
    return
    
def testBlob():
    """write a file as binary data to blob and read back"""
    from ZODB.PersistentMapping import PersistentMapping
    import mimetypes
    from PILView import PILViewer

    DB = PDatabase(server='localhost',port=8090,
              username='farrell',
              password='123')


    def addfile(fname):
        myblob = Blob()
        b=myblob.open('w')
        o=open(fname)
        data = o.read()
        b.write(data)
        print b.name
        b.close()
        return myblob

    dirlist=os.listdir(os.getcwd())

    for f in dirlist:
        m = mimetypes.guess_type(f)[0]
        if m != None and 'image' in m:
            print f
            b=addfile(f)
            DB.add(f)
            DB.data[f]['testfile']=FileRecord(name=f,blob=b)
    DB.commit()
    for k in DB.data:
        if not DB.data[k].has_key('testfile'):
            continue
        rec = DB.data[k]['testfile']
        myblob = rec.blob
        f = myblob.open("r")
        print f.name

    #app = PILViewer(imgfile=f.name)
    #app.mainloop()
    DB.close()
    return

def testRelstorage():
    
    import ZODB, transaction
    from ZODB import FileStorage, DB
    from relstorage.adapters.mysql import MySQLAdapter
    from relstorage.storage import RelStorage
    from MySQLdb import OperationalError

    server='peat.ucd.ie'
    username='guest'
    password='123'
    project='filestest'
    port=8080
    adapter = MySQLAdapter(host=server, user=username,
                                    passwd=password, db=project, port=port)                           
    storage = RelStorage(adapter, shared_blob_dir=False, blob_dir='tempblob')
    db = DB(storage)
    connection = db.open()
    print storage
    connection = db.open()
    dbroot = connection.root()
    data = dbroot['data']
    
    def addfile(fname):
        myblob = Blob()
        b=myblob.open('w')
        o=open(fname)
        data = o.read()
        b.write(data)
        print b.name
        b.close()
        return myblob
    f='gogh.chambre-arles.jpg'
    b=addfile(f)
    data['aaa'] = FileRecord(name=f,blob=b)
    #t = transaction.get()
    #t.commit()
    return

def testLoading():
    """Test loading times for large DB"""   
    from PEATDB.PEATTables import PEATTableModel
    t1 = time()
    DB = PDatabase(local='large.fs')
    #DB = PDatabase(server='localhost', port=8080, username='farrell',
    #                     password='123', project='large')    
    t2=time()
    print round(t2 - t1,2)
    
    print DB.getRecs()
    t3=time()
    print round(t3 - t2,2)
    
    print DB.meta
    t4=time()
    print t4 - t3
    
    M = PEATTableModel(DB)
    print M
    t5=time()
    print t5 - t4
        
    return

def testMemory():
    """test memory behaviour of DB, how can cache be managed"""
    DB = PDatabase(server='localhost', username='farrell',
                         password='123', project='novo')
    #DB = PDatabase(local='large.fs')

    print DB
    db = DB.db
    print db.cacheSize()
    for k in DB.getRecs()[:50]:
        #print k
        r=DB[k]
        r.name
        if db.cacheSize()>500:
            db.cacheMinimize()

    print db.cacheSize()
    return

def importTest():
    DB=PDatabase(local='import.fs')
    DB.importCSV(filename='testdata.csv')

    return

def convertClass():
    """Convert records to proper module name"""

    project='test'
    DB=PDatabase(server='localhost',port=8080,
                 username='farrell',password='123',
                 project=project)

    for k in DB.getRecs():
        r=PEATRecord(k)
        rec = DB.data[k]
        print rec.__class__
        for f in rec.getFields():
            r[f] = rec[f]
        DB.data[k] = r
        print DB.data[k].__class__
    DB.commit(note='convert')
    return

def convertEkinprjs(local=None, server=None,
                    project=None, username=None):
    """convert old ekin prjs in a db to new"""
    if local != None:
        DB = PDatabase(local=local)
    elif server != None:
        DB = PDatabase(server=server, username='farrell', port=8080,
                          password='123', project=project)
    
    for f in DB['userfields']:
        if DB['userfields'][f]['field_type'] in ekintypes:
            print f
            for r in DB.getRecs():
                rec = DB[r]
                if rec.has_key(f):
                    E=rec[f]
                    E.checkDatasets()
                    for d in E.datasets:
                        ek=E.getDataset(d)
                        #ek.prettyPrint()   
                    rec[f] = E                
                  
    print DB.getChanged()  
    DB.commit('converted ekin data')   
    return
    
def remodel():
    #DB=PDatabase(local='hewlsample.fs')
    DB=PDatabase(server='localhost',port=8080,
                 username='farrell',password='123',
                 project='novo')    
    print DB, 'curr ref:', DB.meta.refprotein
    '''rec= DB['wt+D52N']
    print rec.getDisplayAttribute('Structure')
    rec['Structure'] = rec.Structure'''
    
    r='6 c9'
    rec= DB[r]    
    DBActions.checkModels(DB=DB, selected=[r])
    print rec.Structure
    #print rec.aaseq
    '''for r in DB.getRecs():
        print r, type(DB.get(r).Structure)
        print DB.get(r).getDisplayAttribute('Structure')'''

    #DB.commit(note='modelling')
    return

if __name__ == '__main__':
    path=os.getcwd()
    #createdb(local=os.path.join(path,'large.fs'), norecs=5000)
    #createdb(server='localhost', username='zodbuser',
    #                project='large',norecs=10000)
    #importOldProj(datadir='/local/farrell/peat_projects/.TIT_DB.PEAT',
    #                server='localhost', username='farrell', project='titration_db')
    #checkdb(server='peat.ucd.ie', project='titration_db')
    #setDisplayFields(server='localhost', username='farrell',
    #                  project='titration_db')
    #testBlob()
    #speedTest1()
    #print Timeit(speedTest1)
    #testLoading()
    #testnewRecord()
    #testMemory()
    #importTest()
    #remodel()
    #convertEkinprjs(local='tit_db.fs')#server='localhost', username='farrell', project='titration_db')
    testRelstorage()

