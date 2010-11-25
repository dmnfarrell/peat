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

"""Class using ZODB/ZEO as backend"""

import os, sys, types
import numpy
import ZODB, transaction
from ZODB import FileStorage, DB
from BTrees.OOBTree import OOBTree
from ZODB.PersistentMapping import PersistentMapping
from ZODB.blob import Blob
from PEATDB.Record import Record, PEATRecord, Meta, FileRecord
from time import time


class zDatabase(object):
    """Class to replace the peat Database class using the ZODB backend.
       All I/O to disk and server done from this class."""

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
            self.meta = self.dbroot['meta']            
            self.userfields = self.meta.userfields   

            #temp for old dbs - convert meta.info to persistent if it isn't
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
            #self.cellcache = self.dbroot['cellcache']
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
        """Create a BTree object for the root Protein data and a
           seperate object for meta data like userfields.
           All objects added to dbroot should be persistent aware"""
        self.data = self.dbroot['data'] = OOBTree()
        self.meta = self.dbroot['meta'] = Meta('PEAT-DB database')
        self.meta.staticfields = { 'name': 'text', 'Mutations':'text',
                                     'Structure':'PDB'}
        #self.cellcache = self.dbroot['cellcache'] = OOBTree()
        self.blobs = self.dbroot['blobs'] = OOBTree()
        #self.meta.usedisplaycache = False
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

    def setCacheSize(self, s):
        self.meta.info['cachesize'] = s
        return

    def clean(self):
        """Periodic packing to prevent it from consuming your entire disk"""
        self.db.pack()
        return

    def commit(self, user='anon', note='minor'):
        """Commit change in the DB to zodb.
           We treate conflicts as application level errors, so the
           method returns false if it fails to commit"""
        t = transaction.get()
        t.setUser(user)      
        t.note(note)
        from ZODB import POSException
        try:
            t.commit()            
            #self._checkMemCache()
        except POSException.ConflictError:
            print 'there is a conflict'
            print 'following objects changed: ', self.getChanged().keys()            
            return False

        return True

    def abort(self):
        """Abort transactions since last commit"""
        #try:
        transaction.get().abort()
        #except NotImplementedError:
        return

    def supportsUndo(self):
        """Check is storage used supports undo"""
        return self.db.supportsUndo()

    def getData(self):
        return self.data

    def getMeta(self):
        return self.meta

    def recs2List(self):
        return list(self.data.keys())
    
    def listRecs(self):
        for r in self.data.keys():
            print self.data[r]
        return

    def length(self):
        return len(self.data.keys())

    def getSize(self):
        return round(self.db.getSize()/1024.0,1)

    def add(self, key, record=None):
        if key in self.data.keys():
            print 'key %s is already present' %key
            return False
        if record == None:
            record = Record(name=key)
        self._checkMemCache()    
        self.data[key] = record
        self.data._p_changed = 1
        return

    def delete(self, key):
        if not key in self.data:
            print 'no such record'
        else:
            print 'deleted', key            
            del self.data[key]
            self.data._p_changed = 1
        return

    def addField(self, name, fieldtype, action=None):
        """Add a field spec"""
        if name in self.meta.userfields:
            return False
        fieldinfo = {'field_type': fieldtype, 'show': True, 'default_value': ''}
        self.meta.userfields[name] = fieldinfo
        self.meta._p_changed = 1
        return True

    def deleteField(self, fname, callback=None):
        """Remove the field from all records"""
        if not fname in self.meta.userfields:
            print 'field name not in userfields..'
        i=0    
        for k in self.data.keys():
            #we cannot minimize cache between transactions, so 
            #may need to manually commit in batches to prevent filling memory
            '''if self.db.cacheSize() > self.meta.info['cachesize']:
                self.commit('field deletion')
                self.db.cacheMinimize()'''
            print self.db.cacheSize() 
            if fname in self.data[k].keys() and type(self.data[k]) is PEATRecord:                
                #del self.data[k][fname]
                self.data[k].delete(fname)
            i+=1    
            if callback != None:               
                callback(float(i)/self.length()*100)
        del self.meta.userfields[fname]
        self.meta._p_changed = 1
        return

    def showField(self, name, status=True):
        """Hide from table"""
        if name in self.meta.userfields:
            self.meta.userfields[name]['show'] = status
            self.meta._p_changed = 1
        return

    def getFields(self):
        """Get fields"""
        return self.meta.userfields.keys()

    def getSimpleFields(self):
        """Get fields with only text or number data"""
        fields = []
        for f in self.getFields():
            if self.getFieldType(f) == 'text':
                fields.append(f)            
        return fields
    
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
        self._checkMemCache()
        return self.data[key]

    def getRecs(self):
        return self.data.keys()

    def getRecordName(self, name):
        """Get actual record name from name field (might be different)"""
        for r in self.getRecs():
            if name == self[r].name:
                return r
        return None
    
    def __getitem__(self, key):
        """Allow metafields to be retrieved from DB as keys"""
        self._checkMemCache()
        if key in Meta.special:
            return self.meta[key]
        return self.data[key]

    def __repr__(self):
        return 'DB with %s records' %len(self.data.keys())

    def importDict(self, importdata, namefield='name', overwrite=True):
        """Import list of dicts, each one is a record"""
        if len(importdata) == 0:
            return        
        fields = importdata[0].keys()       
        if not namefield in fields:
            print 'no such field for keyname field'
            namefield = fields[0]
            #return
        if namefield == 'name': fields.remove(namefield)
        for f in fields:
            if f not in self.meta.staticfields:
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
        persisted = [self.meta, self.data, self.meta.info]
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

class PDatabase(zDatabase):
    """Class to replace the peat Database class. Inherits from zDatabase.
       with protein stuff added and has a representation in PEATTableModel"""

    def __init__(self, server=None, local=None, storage=None, backend='mysql',
                 port=3306, project='test',
                 username=None, password=None, blob_dir=None):
        zDatabase.__init__(self, server, local, storage, backend,
                            port, project, username, password, blob_dir)
        self.ekintypes = ['General', 'NMR titration', 'Protein Stability',
                      'pH Activity profile', 'Michaelis-Menten kinetics',
                      'pH stability profile', 'Amyloid formation']
        self.refprotein = None
        return

    def add(self, key, record=None):
        if key in self.data.keys():
            print 'key %s is already present' %key
            return False
        if record == None:
            record = PEATRecord(name=key)
        self._checkMemCache()    
        self.data[key] = record
        self.data._p_changed = 1
        return
    
    def getRecordWithMutation(self, mutationstring):
        """Get first record with the provided mutation"""
        for p in self.getRecs():
            prot=self.DB.get(p)
            if mutationstring == prot.Mutations:
                return p
        return None
    
    def getTextFields(self):
        flds = ['name']
        for f in self.getFields():
            if self.meta.userfields[f]['field_type'] == 'text':
                flds.append(f)
        return flds
    
    def getDictFields(self):
        flds = []
        #get fields that we store as dicts
        dicttypes = ['dict', 'Link', 'Table', 'Sequence', 'Notes']
        for f in self.getFields():
            if self.meta.userfields[f]['field_type'] in dicttypes:
                flds.append(f)
        return flds
    
    def getEkinFields(self):
        flds = []
        for f in self.getFields():
            if self.meta.userfields[f]['field_type']  in self.ekintypes:
                flds.append(f)
        return flds
        
    def addLabbook(self, name, data=None):
        """Create an extra labbook
           We store these in meta.labbooks"""
        if not hasattr(self.meta, 'labbooks'):
            self.meta.labbooks = PersistentMapping()
        self.meta.labbooks[name] = data
        return
    
    def createLabbookSheet(self, sheet, model=None):
        if model == None:
            from TableModels import TableModel
            model = TableModel()
        self.meta.labbook[sheet] = model.getData()
        return model
          
    def getLabbookSheet(self, sheet):
        """Get a Labbook sheet model instance"""
        data = self.meta.labbook[sheet]
        from TableModels import TableModel
        model = TableModel(data)
        return model

    def getLabbookData(self, sheet, fields):
        """Get Labbook column data"""
        S = self.meta.labbook[sheet]
        from TableModels import TableModel
        model = TableModel(S)
        data={}
        for f in fields:
            data[f] = model.getColumnData(columnName=f)
        return data

    def saveLabbook(self, sheet, model):
        """Save a model instance into the provided sheet"""
        tmp = self.meta.labbook
        tmp[sheet] = model.getData()
        self.meta.labbook = tmp
        return
       
    def saveLabbooktoFile(self, filename):
        """Save labbook to a file"""
       
        fd=open(filename,'w')
        import pickle
        pickle.dump(self.meta.labbook,fd)
        fd.close()        
        return        
    
    def addDNAseq(self, protein_name, DNAseq):
        """Add a DNA sequence to a protein"""
        if self.isknownProtein(protein_name):
            # Do we already have a DNA sequence?
            self.data[protein_name].DNAseq = DNAseq
            return True
        return None

    def addProtseq(self,protein_name,sequence,start_aanumber=1,
                        seq_format='detect',update=True):
        """
         Add the protein sequence - blind faith for now...
         We add this protein sequence as currentres. If the protein has a parent, then it must be set
         elsewhere
        """
        if seq_format=='detect':
            import types
            if type(sequence)==types.ListType:
                if len(sequence[0])==3:
                    pass
                elif len(sequence[0])==1:
                    sequence=self.convert_list1_2_list3(sequence)
            elif type(sequence)==types.StringType:
                status=self.convert_string3_2_list3(sequence)
                if status:
                    sequence=status
                else:
                    sequence=self.convert_string1_2_list3(sequence)
        else:
            actions={'1AA_list':self.convert_list1_2_list3,
                     '3AA_list':None,
                     '1AA_string':self.convert_string1_2_list3,
                     '3AA_string':self.convert_string3_2_list3}
            if actions.has_key(seq_format):
                action=actions[seq_format]
                if action:
                    sequence=action(sequence)
                else:
                    print 'Unknown sequence format'
                    raise Exception()

        # We have the sequence in the right format, so simply insert it

        if self.isknownProtein(protein_name):
                # Create the sequence and add it to the DB
                count=start_aanumber
                seq=[]
                import string
                for aa in sequence:
                    number=':'+string.zfill(count,4)
                    seq.append([number,aa])
                    count=count+1
                self.data[protein_name].aaseq = seq

                # Add the start aa number
                self.data[protein_name].aa_number_start = start_aanumber
                return True
        else:
            return None

    def isknownProtein(self,name):
        """Does the database already contain a protein with this name?"""
        if self.data.has_key(name):
            if self.data[name].type=='protein' or self.DB[name].type=='enzyme':
                return True
        return False

    def addProperty(self,protein_name,category,this_dict):
        """Add the properties in dict to the information for protein_name.
        category is the key/name that the info will be associated with"""

        if self.isknownProtein(protein_name):
            if type(this_dict) is type({}):
                self.data[protein_name][category]=this_dict.copy()
            else:
                self.data[protein_name][category]=this_dict
            
            return 1
        return None

    def convert_string1_2_list3(self,sequence):
        """Convert a string of AA1 to a list of AA3"""
        newseq=[]
        import DNAtool.mutation
        for letter in sequence:
            newseq.append(DNAtool.mutation.one_to_three[letter])
        return newseq

    def convert_string3_2_list3(self,sequence):
        """Convert a string of AA3 to a list of AA3"""
        newseq=[]
        thisaa=''
        import DNAtool.mutation
        for letter in sequence:
            thisaa=thisaa+letter
            if len(thisaa)==3:
                if DNAtool.mutation.three_to_one.has_key(thisaa):
                    newseq.append(thisaa)
                else:
                    return None
                thisaa=''
        return newseq

    def get_next_freename_newprot(self):
        """Get the next free name in the database"""
        import string
        found=None
        base='E'
        count=1
        while not found:
            name='%s%s' %(base,string.zfill(count,3))
            if not self.data.has_key(name):
                found=name
            count=count+1
        return found

    #this method seems redundant to me..
    def newmolecule(self,protein_name,update=True):
        """
        This function files a new molecule in the Database
        protein_name is the name of the molecule (DB identifier)
        Enter all the data in dict into the database
        """
        if self.isknownProtein(protein_name):
            return None
        else:
            #self.DB[protein_name]=PEAT_dict(changelist=self.changelist,protein=protein_name)
            self.add(protein_name)

        return True

    def get_AA_sequence(self,protein_name):
        return self.get_protein_sequence(protein_name)

    def get_protein_sequence(self, protein_name):
        """Get the amino acid sequence for this protein"""
        if self.data[protein_name].has_key('aaseq'):
            aaseq=self.data[protein_name]['aaseq']
        else:
            aaseq=None
        return aaseq

    #these should be at protein rec level?

    def storePDB(self, name, X, AlignmentMap=None):
        """Store the PDB file for the selected protein"""
        rec = self.data[name]
        rec.Structure = X.writepdb('dummy',nowrite=1)
        rec.Structuretype = 'xray'
        # Store the map between sequences
        if AlignmentMap != None:
            rec.Structure_alnseq_EATrecord = AlignmentMap['OrigAa']
            rec.Structure_alnseq_PDBfile = AlignmentMap['AlignedAa']
        return

    def getStructure(self, protein, field_name=None):
        """Get a structure for this record.
           Return as a list of pdblines,Protool instance"""
        if not self.data.has_key(protein):
            return None,None
        if field_name == None:
            field_name = 'Structure'
        rec = self.get(protein)
        S = rec.Structure
        if type(S) is types.ListType:
            # Real X-ray structure
            import Protool
            X = Protool.structureIO()
            pdblines = S
            X.parsepdb(pdblines)
            return pdblines, X
        
        elif type(S) is type({}):
            # Otherwise model on the fly
            if self.data[protein][field_name].has_key('Rotamer_operations'):
                from PE import ProteinEngineering
                pdblines, X = ProteinEngineering.modelontheFly(self, S)
                if not pdblines:
                    import tkMessageBox
                    tkMessageBox.showwarning('Could not model structure',
                                                 'Something went wrong during the modelling process.')
                    return None,None
                return pdblines,X
        return None,None


