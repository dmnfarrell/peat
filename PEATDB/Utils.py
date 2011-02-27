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

"""Factory methods for handling PEAT Datbases, such as combining
   two databases"""

def loadDB(prj, remote=True, settings={}):
    """Load a local or remote db, settings are the remote db arguments
       such as server, username, password and project name"""
    if remote==True:
        DB = PDatabase(**settings)
    else:
        DB = PDatabase(local=prj)
    return DB

def combineDBs(db1, db2):
    """Combine tow databases"""
    comb = PDatabase()
    return comb

def createDBonServer(self, settings={}, callback=None):
    """Create a project on the mysql server"""
    import MySQLdb as mysql
    db = mysql.connect(**settings)
    c = db.cursor()
    cmd = "create database " + dbname + ";"
    c.execute(cmd)
    return True
    
def copyDBtoServer(localdb, remotedb):
    """Copy the contents of a local db to a remote one, requires
      that user has already made the remote connection and passes
      the remote db object"""
    
    return True
                   
