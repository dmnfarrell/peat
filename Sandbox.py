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

"""test ZODB stuff"""

import os, sys, types
import numpy
import ZODB, transaction
from ZODB import FileStorage, DB
from BTrees.OOBTree import OOBTree
from persistent import Persistent
from ZODB.PersistentMapping import PersistentMapping
from ZODB.PersistentList import PersistentList
from ZODB.blob import Blob
from PEATDB.Record import PEATRecord, Meta, FileRecord
from time import time

class Record1(object):
    """Holds the field info and some other attributes for a record
       Dynamic reference to the data in the db, so we can handle
       records without an entire 'record' object"""

    def __init__(self, db, name, fields):
        self.__dict__['db'] = db
        self.__dict__['name']=name
        self.__dict__['fields']=fields        
        return

    def fields(self):
        return self.fields
    
    def __getitem__(self, key):
        if key == 'name' or key == 'Name':
            return self.name
        return self.db.data[self.name+'_'+key]

    def __setitem__(self, key, data):
        if not key in self.fields:
            self.fields.append(key)
        self.db.data[self.name+'_'+key] = data       
        return

    def __getattr__(self, item):
        """Maps values to attributes.
        Only called if there *isn't* an attribute with this name
        """
        try:
            return self.__getitem__(item)
        except KeyError:
            raise AttributeError(item)

    def __setattr__(self, key, data):
        """Maps attributes to key values."""      
        return self.__setitem__(key, data)

    def __str__(self):
        return 'Record with %s fields' %len(self.fields)
    
    
class myDatabase(object):
    """New version of peat db that distributes objects by rec/field item
       rather than just per each record.
       So we store each cell thereby rarely need to touch all objects in the DB.
       """

    backends = ['mysql', 'zeo']

    def __init__(self, server=None, local=None, storage=None, backend='mysql',
                 port=3306, project='test',
                 username=None, password=None, blob_dir=None):
        """Start with local or server connection to zodb/zeo"""
        import logging
        logging.basicConfig()
        self.project = project
        self.backend = backend
        self.errormessage = None
        if blob_dir==None:
            self.setBlobdir(relstorage=True)
        else:
            self.blob_dir = blob_dir

        if storage != None:
            self.storage = storage
        elif server != None:
            result = self.createStorage(server, username, password,
                                         project, port, backend)
            self.storagetype = 'remote'
            if result == False:
                return None
        elif local != None:
            #get the blob dir for the local filename
            localblobdir = self.getLocalBlobdir(local)
            self.storage = FileStorage.FileStorage(local,
                                                   blob_dir=localblobdir)          
            self.storagetype = 'local'
        else:
            return None
        self.db = DB(self.storage)
        self.db.setCacheSize(50)
        self.connection = self.db.open()
        self.dbroot = self.connection.root()
        self.open = True

        #define data
        if self.dbroot.has_key('data'):
            self.data = self.dbroot['data']
            self.records = self.dbroot['records']
            self.meta = self.dbroot['meta']           
            self.userfields = self.meta.userfields
           
            if not self.meta.has_key('info'):
                self.meta.info = PersistentMapping({'project': '', 'users':'', 'cachesize':100})
                self.commit()
            if type(self.meta.info) != PersistentMapping:
                x = PersistentMapping(self.meta.info)
                self.meta.info = x
                self.commit(note='convert db info')
            if not self.meta.info.has_key('cachesize'):
                self.setCacheSize(200)
                self.commit(note='cachesize set')
            self.cellcache = self.dbroot['cellcache']
            try:
                self.blobs = self.dbroot['blobs']
            except:
                self.blobs = self.dbroot['blobs'] = OOBTree()
                self.commit(note='blob dir')
        else:
            self.createData()

        return

    def createStorage(self, server, username, password,
                      project, port, backend):
        """Create the specified storage from given options"""

        if backend == 'mysql':
            #we now use RelStorage by default, mysql backend
            from relstorage.adapters.mysql import MySQLAdapter
            from relstorage.storage import RelStorage
            from MySQLdb import OperationalError
            try:
                self.adapter = MySQLAdapter(host=server, user=username,
                                            passwd=password, db=project, port=port,
                                            connect_timeout=3)
                self.storage = RelStorage(self.adapter)
            except OperationalError, message:
                self.connection = None; self.dbroot = None
                self.errormessage = message
                return False
        else:
            #otherwise we try to use clientstorage/ZEO
            from ZEO.ClientStorage import ClientStorage
            server_and_port = (server, port)
            self.storage = ClientStorage(server_and_port,
                                             storage=project,
                                             username=username,
                                             password=password,
                                             wait=False,
                                             wait_timeout=50,
                                             blob_dir=self.blob_dir)
            #can't raise exception here if connection fails..
            if 'disconnected' in self.storage.getName():
                self.connection = None; self.dbroot = None
                self.errormessage = 'Zeo error'
                return False
        return True

    def setBlobdir(self, blob_dir=None, relstorage=True):
        """Set the blob dir for use with clientstorage"""
        home = os.path.expanduser("~")
        if relstorage == True:
            self.blob_dir = os.path.join(home, '.peatblob')
        else:
            self.blob_dir = os.path.join(home, '.zblob')
        return

    def getLocalBlobdir(self, path):
        """One per local file"""
        d = os.path.abspath(path)+'.blob'
        return d

    def createData(self):
        """Create a OOBTree object for the root Protein data and a
           seperate object for meta data like userfields.
           All objects directly added to dbroot should be persistent aware
           New db adds each cell individually, not whole rec structure
        """
        self.data = self.dbroot['data'] = OOBTree()
        #self.records = self.dbroot['records'] = PersistentList()
        self.records = self.dbroot['records'] = OOBTree()
        #meta is used for all other non-record stuff
        self.meta = self.dbroot['meta'] = Meta('PEAT-DB database')
        self.meta.staticfields = { 'Name': 'text', 'Mutations':'text',
                                     'Structure':'PDB'}
        self.cellcache = self.dbroot['cellcache'] = OOBTree()
        self.blobs = self.dbroot['blobs'] = OOBTree()
        self.meta.usedisplaycache = False
        #use meta.info for any misc settings specific to the DB
        self.meta.info = PersistentMapping({'project': '', 'users':'', 'cachesize':100})
        #try to limit no of objects to store in RAM
        self.meta.cachesize = 100
        self.commit()
        return

    def close(self, commit=False):
        """Close connection, puts DB offline"""

        if commit == True:
            print 'commiting'
            transaction.get().commit()
        else:
            transaction.get().abort()

        #self.dbroot.clear()
        self.connection.close()
        self.db.close()
        self.storage.close()
        self.open = False
        return

    def createCellCache(self):
        """Create a cache of the displayed values for each record/field
           Primarily used for quick display in the table"""
        print '(re)creating db display cache'
        fields = self.meta.userfields.keys() + self.meta.staticfields.keys()
        self.cellcache.clear()
        for f in fields:
            self.cellcache[f] = {}
            if f in self.meta.userfields:
                ftype = self.meta.userfields[f]['field_type']
            else:
                ftype = self.meta.staticfields[f]
            for r in self.records.keys():
                p = self.data[r]
                self.cellcache[f][r] = p.getDisplayAttribute(f, ftype)
        self.cellcache._p_changed = 1
        self.meta.usedisplaycache = True
        return

    def clearCellCache(self):
        """Clear display cache """
        print 'reset db display cache'
        self.cellcache.clear()
        self.meta.usedisplaycache = False
        return

    def updateCellCache(self, name, field):
        """Update cache for cell display whenever we add/remove/change fields"""
        if field in self.meta.userfields:
            ftype = self.meta.userfields[field]['field_type']
        else:
            ftype = self.meta.staticfields[field]
        if not self.cellcache.has_key(field):
            self.cellcache[field] = {}
        p = self.data[name]
        self.cellcache[field][name] = p.getDisplayAttribute(field, ftype)
        self.cellcache._p_changed = 1
        return

    def setCacheSize(self, s):
        self.meta.info['cachesize'] = s
        return

    def clean(self):
        """Periodic packing to prevent it from consuming your entire disk"""
        self.db.pack()
        return

    def commit(self, user='anon', note=None):
        """Commit change in the DB to zodb.
           We treate conflicts as application level errors, so the
           method returns false if it fails to commit"""
        t = transaction.get()
        t.setUser(user)
        if note==None:
            note='minor'
        t.note(note)
        from ZODB import POSException
        try:
            t.commit()
        except POSException.ConflictError:
            print 'there is a conflict'
            print 'following objects changed: ', self.getChanged().keys()
            #t.abort()
            return False

        return True

    def abort(self):
        """Abort transactions since last commit"""
        #try:
        transaction.get().abort()
        #except NotImplementedError:
        return

    def supportsUndo(self):
        """Check if storage used supports undo"""
        return self.db.supportsUndo()

    def getData(self):
        return self.data

    def getMeta(self):
        return self.meta

    def listRecs(self):
        for r in self.records.keys():
            print self.records[r]
        return

    def length(self):
        return len(self.records.keys())

    def getSize(self):
        return round(self.db.getSize()/1024.0,1)

    def add(self, key):
        if key in self.records:
            print 'key %s is already present' %key
            return False
        self.records[key] = []
        #self.records.append(key)
        self._checkMemCache()        
        self.records._p_changed = 1
        return

    def delete(self, key):
        if not key in self.records:
            print 'no such record'
        else:
            print 'deleted', key            
            del self.records[key]
            #remove from data also
            
            self.records._p_changed = 1
            self.data._p_changed = 1
        return

    def addField(self, name, fieldtype, action=None):
        """Add a field spec"""
        if name in self.meta.userfields:
            return False
        fieldinfo = {'field_type':fieldtype, 'action':action, 'default_value': ''}
        self.meta.userfields[name] = fieldinfo
        self.meta._p_changed = 1
        return True

    def deleteField(self, fname, callback=None):
        """Remove the field from all records"""
        if not fname in self.meta.userfields:
            print 'field name not in userfields..'
        i=0    
        for k in self.data.keys():
            self._checkMemCache()
            print self.db.cacheSize() 
            if fname in self.data[k].keys() and type(self.data[k]) is PEATRecord:
                self.data[k].delete(fname)
            i+=1    
            if callback != None:               
                callback(float(i)/self.length()*100)
        del self.meta.userfields[fname]
        self.meta._p_changed = 1
        return

    def hideField(self, name):
        """Remove from userfields, but keep in records"""
        if name in self.meta.userfields:
            del self.meta.userfields[name]
        return

    def getFields(self):
        """Get fields"""
        return self.meta.userfields.keys()

    def getFieldType(self, field):
        """Get field type"""
        if self.meta.userfields.has_key(field):
            return self.meta.userfields[field]['field_type']

    def addBlob(self, recname, field, filename):
        """Adds an external file as a blob to the DB"""
        f=os.path.split(filename)[1]
        self.add(recname)
        blob = self.getBlob(filename)
        self.data[recname][field] = FileRecord(name=filename, blob=blob)
        return

    def getBlob(self, filename):
        """Create a blob from a file"""
        from ZODB.blob import Blob
        myblob = Blob()
        b=myblob.open('w')
        o=open(filename)
        data = o.read()
        b.write(data)
        b.close()
        return myblob

    def undo(self, id, user, note, callback=None):
        """Undo a commit, given id"""
        from ZODB import POSException
        try:
            self.db.undo(id=id)
            self.commit(user=user, note=note)
        except POSException.UndoError, POSException.MultipleUndoErrors:
            print 'failed, later commits changed the same data'
            callback()
            self.abort()
        return

    def undoLog(self):
        return self.storage.undoLog(0, sys.maxint)

    def pack(self,days=7):
        """Pack DB"""
        if self.storagetype == 'remote' and self.backend == 'mysql':
            #pack relstorage
            import ZODB.serialize, time
            t = time.time() - float(days) * 86400.0
            self.storage.pack(t, ZODB.serialize.referencesf)
        else:                
            self.db.pack(days=days)
        return

    def get(self, key):
        return self[key]

    def getRecs(self):
        #return list(self.records.keys())
        return self.records
    
    def __getitem__(self, key):
        """Allow metafields to be retrieved from DB as keys"""
        self._checkMemCache()
        if key in Meta.special:
            return self.meta[key]
        if key in self.records:
            fields = self.records[key]
            rec = Record1(self, key, fields)
         
        return rec

    def __repr__(self):
        #return 'DB with %s records' %len(self.records.keys())
        return 'DB with %s records' %len(self.records)
    
    def importDict(self, importdata, namefield='name', overwrite=True):
        """Import list of dicts, each one is a record"""
        if len(importdata) == 0:
            return
        
        fields = importdata[0].keys()
        if not namefield in fields or not 'name' in fields:
            print 'no such field for keyname field'
            namefield = fields[0]
            #return
        fields.remove(namefield)
        for f in fields:
            self.addField(f, 'text')
        for d in importdata:
            name = d[namefield]
            self.add(name)
            for f in fields:
                self.data[name][f] = d[f]        
        #print self
        return

    def importCSV(self, filename, namefield='name'):
        """Import from a CSV file"""
        from IO import Importer
        IM = Importer()
        importdata = IM.doImport(filename)
        self.importDict(importdata, namefield)
        return

    def addFile(self, blob):
        self.connection.add(blob)
        return

    def isChanged(self):
        """Check if DB has been changed"""
        persisted = [self.meta, self.data, self.cellcache, self.meta.info]
        for p in persisted:
            if p._p_state == 1:
                return True
        for k in self.meta.keys():                    
            if not type(self.meta[k]) is PersistentMapping:
                continue
            if self.meta[k]._p_state == 1:
                return True
        for k in self.data.keys():
            if self.data[k]._p_state == 1:
                return True
        return False

    def getChanged(self):
        """Get dict of changed objects"""
        changed = {}
        if self.meta._p_state == 1:
            changed['meta'] = self.meta
        for k in self.meta.keys():                    
            if not type(self.meta[k]) is PersistentMapping:
                continue
            if self.meta[k]._p_state == 1:                
                changed[k] = 1
        if self.cellcache._p_state == 1:
            changed['cellcache'] = self.cellcache
        if self.data._p_state == 1:
            changed['data'] = True  #we can't cache the entire data tree...
        for k in self.data.keys():
            if self.data[k]._p_state == 1:
                changed[k] = self.data[k]
        return changed

    def getExtFilenames(self):
        """Return a list of all external files stored in DB"""
        files = {}
        for field in self.getFields():
            if self.getFieldType(field) == 'File':
                for k in self.getRecs():
                    if self.get(k).has_key(field):
                        rec = self.get(k)[field]
                        f = rec.getFile()
                        name = os.path.basename(rec.filename)
                        files[name] = f.name
        return files

    def _checkMemCache(self):
        """Need to check the memory isn't being overloaded
           Call this when addressing lots of recs serially, eg. inside getitem
           not very fine grained"""
        if self.db.cacheSize() > self.meta.info['cachesize']:
            self.db.cacheMinimize()
        return

def test():
    """Test IO and mem usage"""
    import copy
    DB = myDatabase(local='newdb.fs')
    #DB = PDatabase(local='olddb.fs')
    db = DB.db
    from PEATDB.Ekin.Base import EkinProject
    import copy
    DB.addField('stab', 'text')
    DB.addField('ekin', 'General')
    DB.addField('ekin2', 'General')
    DB.meta.info['cachesize'] = 50
    E = EkinProject()
    E.openProject('big')  

    '''for i in range(200):
        Ec =  copy.deepcopy(E)           
        name = 'rec'+str(i)        
        DB.add(name)        
        r = DB[name]
        r['ekin'] = Ec
        r.ekin2 = i
        print r
    print DB   
    DB.commit()'''
   
    #DB.deleteField('ekin')
    
    '''for i in DB.getRecs():
        #x=DB[i]['ekin']
        #y=DB[i]['ekin2']
        print i ''' 
    print DB, list(DB.data.keys())  
    #print db.cacheSize()
 
    return
    
if __name__ == '__main__':
    path=os.getcwd()
    test()
