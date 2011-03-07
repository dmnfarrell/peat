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
# Author: Damien Farrell 2011

import pickle, sys, os, copy, time, types
import numpy
from PEATDB.Base import PDatabase 

"""Helper methods for handling PEAT Datbases, such as copying, merging
   two databases"""

settings={'server':'localhost','username':'guest',
           'password':'123'}

def loadDB(prj, remote=False, settings={}):
    """Load a local or remote db, settings are the remote db arguments
       such as server, username, password and project name"""
    if remote==True:
        DB = PDatabase(project=prj,**settings)
    else:
        DB = PDatabase(local=prj)
    return DB

def mergeDBs(DB1, DB2):
    """Combine two databases"""
    newDB = PDatabase()
    return newDB

def copyDB(DB1, DB2, overwrite=True):
    """Copy one db to another"""
    count=0
    for r in DB1.getRecs():
        rec = copy.deepcopy(DB1.data[r])
        DB2.add(r, record=rec)
        count+=1
    for m in DB1.meta.keys():
        DB2.meta[m] = copy.deepcopy(DB1.meta[m])
    print 'copied db successfully'
    return DB2

def saveDBCopy(DB, filename, callback=None):
    """Save local copy of a remote or another local DB"""
    import copy
    total = len(DB.getRecs())
    if filename == '' or filename == None:
        return False
    if os.path.exists(filename):
        for i in ['.lock','.index','']:
            os.remove(filename+i)
    newDB = PDatabase(local=filename)
    newDB = copyDB(DB, newDB)
    newDB.commit()    
    newDB.close() 
    return newDB

def createDBonServer(prj, settings={}, callback=None):
    """Create a project on the mysql server"""
    import MySQLdb as mysql
    db = mysql.connect(**settings)
    c = db.cursor()
    cmd = "create database " + prj + ";"
    c.execute(cmd)
    return True
    
def copyDBtoServer(DB, prj, settings):
    """Copy the contents of a local db to a remote one, requires
      settings to make the remote connection and get
      the remote db object"""
    remoteDB = loadDB(prj, remote=True, settings=settings)
    print 'source db: %s, remote db: %s' %(DB, remoteDB)
    newDB = copyDB(DB, remoteDB)
    newDB.commit(note='copy')
    return True
                   
if __name__ == '__main__':
    #test    
    db = loadDB(os.path.join('scripts','testing.fs'))
    print db
    #db1 = saveDBCopy(db, filename='test1')
    copyDBtoServer(db, 'test1', settings)
    
