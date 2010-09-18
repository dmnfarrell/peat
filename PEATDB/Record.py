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

import ZODB
from persistent import Persistent
from BTrees.OOBTree import OOBTree
import os, types
from PEATDB.Ekin.Base import EkinProject

class Item(Persistent):
    """Items inside a record have to be put into this class so we can
       have lazy loading of keys in a record"""
    def __init__(self, data):
        self.data = data      
        return
        
class Record(OOBTree):
    """Proto class to represent a persistent object in PEAT.
       We subclass this for protein records and meta info"""
    def __init__(self, name):
        OOBTree.__init__(self)      
        self._display = OOBTree()
        self.name=name
        return
            
    def setDisplayAttribute(self, key, data):
        """Set the text to display for this attribute of the record
           Override this for specific field objects"""
        text = ''
        if type(data) is types.StringType:
            text = data
        else:
            text = ''
        if text == None:
            return ''
        return text

    def getDisplayAttribute(self, key):
        """Get the text to display for this attribute of the record
           This is used for the table display to cache the values"""
       
        if not self._display.has_key(key):
            return ''       
        return self._display[key]
        
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
        if key == 'display':
            return OOBTree.__getitem__(self, key)
        else:
            return OOBTree.__getitem__(self, key).data  

    def __setitem__(self, key, data):  
        if key == 'display':
            OOBTree.__setitem__(self, key, data)
        else:    
            text = self.setDisplayAttribute(key, data)            
            item = Item(data)  
            OOBTree.__setitem__(self, key, item)
            self._display[key] = text
        self._p_changed = 1    
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
        
        
class PEATRecord(Record):
    """We can add protein specific stuff here """

    ekinfields = EkinProject.ekinfields + ['mode']
    ekintypes = ['General', 'NMR titration', 'Protein Stability',
                      'pH Activity profile', 'Michaelis-Menten kinetics',
                      'pH stability profile', 'Amyloid formation']

    def __init__(self, name):        
        Record.__init__(self, name)
        kys = ['Structure','Mutations','aaseq','DNAseq','type']
        for k in kys:
            self[k] = None
        self.type = 'protein'
        return

    def isEkin(self, key):
        if not self.has_key(key):
            return None
        if type(self[key]) is EkinProject:
            return True
        else:
            return False

    def hasStructure(self):
        """Determine if this record has a structure"""
        import types
        if self.Structure == None:
            return 'not available'
        if self.has_key('structuretype'):
            return self.structuretype
        #we need to get rid of this guessing bit.
        if type(self.Structure) is types.ListType and len(self.Structure)>25:
            return 'available'
        elif type(self.Structure) is types.StringType:
            return self.Structure
        elif type(self.Structure) is types.DictType:
            if self.Structure is types.DictType:
                return 'PEAT Model'
            else:
                return 'not available'
        else:
                return 'not available'
        return None    

    def checkStructure(self, structure):
        import types      
        if structure == None:
            return 'not available'     
        elif self.has_key('structuretype'):
            return self.structuretype
        elif type(structure) is types.ListType and len(structure)>25:
            return 'available'
        elif type(structure) is types.StringType:
            return structure
        elif type(structure) is types.DictType:           
            return 'PEAT Model'
        else:
            return 'not available'         
    
    def getMutationString(self):
        """Get reduced mutation string from mutationset if available"""
        mut = self.Mutations
        if mut != None:
            return '+'.join(mut.mutationCodes(reduced=True))
        else:
            return None   
    
    def getAncestry(self, parent, immediate=None):
        """Find out if a record is parent of another record using
           sequence information """
        import sequence_alignment
        from Sequence import SequenceOperations as SQ
        full_record_sequence = self.aaseq       
        full_parent_sequence = parent.aaseq
                
        if parent.has_key('Structure_alnseq_PDBfile') and parent.has_key('Structure_alnseq_EATrecord'):
            parent_aln_PDBfile=parent.Structure_alnseq_PDBfile + '*'
            parent_aln_record=parent.Structure_alnseq_EATrecord + '*'
        else:
            parent_aln_PDBfile=record_sequence,ignored_res=sequence_alignment.Protool2pir(full_record_sequence)
            parent_aln_record=parent_aln_PDBfile

        if full_record_sequence and full_parent_sequence:            
            record_sequence,ignored_res = sequence_alignment.Protool2pir(full_record_sequence)
            parent_sequence,ignored_parent = sequence_alignment.Protool2pir(full_parent_sequence)
            
            # First try the simple option
            
            operations = [None]
            if len(record_sequence) == len(parent_sequence):
                operations = SQ.findSequenceDifferences(record_sequence, parent_sequence,
                                                        full_parent_sequence,
                                                        PDBaln=parent_aln_PDBfile,
                                                        recordALN=parent_aln_record)
                
            if len(operations)>10 or len(record_sequence) != len(parent_sequence):
                NW_align = sequence_alignment.NW(record_sequence,parent_sequence)
                al_seq, al_parent, map_seq, map_parent = NW_align.Align(verbose=True)
                operations = SQ.findSequenceDifferences(al_seq,al_parent, full_parent_sequence,
                                                       PDBaln=parent_aln_PDBfile,
                                                       recordALN=parent_aln_record)
            if operations == [None]:
                raise Exception()
  
            if len(operations) >= 0 and not immediate:
                return True, operations
            elif len(operations) > 0 and immediate:
                if len(operations) == 1:
                    return False, operations
            else:
                return False, operations
        return False, False

        
    def setDisplayAttribute(self, key, data):
        """Set the text to display for this attribute of the record
           This is used for the table display to cache the values"""
        text = ''      
        if key == 'Structure':           
            return self.checkStructure(data)
        if data == None:
            return ''
        if key == 'ORF_selected':
            return ''
        if type(data) is types.StringType:
            text = data
        elif type(data) is types.FloatType or type(data) is types.IntType:
            text = str(data)
        elif type(data) is EkinProject:
            #numkeys = data.length
            numkeys = len(data.datasets)
            text = str(numkeys) + ' recs'         
        elif type(data) is types.DictType:
             if data =='':
                 text = '-'
             elif data.has_key('name'):
                 text = data.name
             elif data.has_key('link'):
                 return data                 
             elif data.has_key('text'):
                 text = ''
                 t = data['text']
                 text = t[0:25]
                 text = text.replace ( "\n", " " )
                 text = text.strip()            
             elif data[data.keys()[0]].has_key('columnnames'):
                 numkeys = len(data)
                 text = str(numkeys) + ' sheets'
             else:
                 text = ''
               
        return text

    def __repr__(self):
        """Return a string representation of all data"""
        s = 'Record name %s with %s fields' %(self.name,len(self.keys()))
        fieldnames = ''
        for key in self.keys(): fieldnames += key+' '
        s += '\nFields:%s' % fieldnames
        return s

class Meta(OOBTree):
    """Class to represent meta info fields in PEAT"""
    special = ['userfields','staticfields','labbook','table','DNAtool_primers']
    def __init__(self, name=None):        
        OOBTree.__init__(self) 
        self.name = name
        for k in self.special:
            self[k] = {}
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
        
    def __repr__(self):
        """Return a string representation of all data"""
        s = 'Meta Info: '
        for key in self.keys():
            s += key+','
        return s


class FileRecord(Meta):
    """Class to store files, these contain blobs with some meta data"""
    def __init__(self, name, blob=None):
        import mimetypes
        self.filename = name
        self.name = os.path.basename(name)
        self.ext = os.path.splitext(name)[1]
        self.size = os.path.getsize(name)
        mimetypes.init()
        self.mimetype = mimetypes.guess_type(name)[0]
        self.blob = blob
        return

    def getFile(self):
        """Get the binary file kept in the blob"""
        if self.blob == None:
            print 'no file'
            return
        f = self.blob.open("r")
        return f

    def __repr__(self):
        """Return a string representation of all data"""
        return str('name:'+self.name+' type:'+self.mimetype)

