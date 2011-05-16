#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-
# All rights reserved
#
# Written by Damien Farrell, Feb 2010

import __main__
__main__.pymol_argv = [ 'pymol', '-qc']
import pymol
import pymol.cmd as cmd
from pymol import stored

from PEATDB.Ekin.Titration import TitrationAnalyser
from PEATDB.Base import PDatabase
from PEATDB.Ekin.Base import EkinProject
from PEATDB.DictEdit import DictEditor
import os
import pickle

cols =  ['15N NMR', '1H NMR']
nuclnames = {'1H NMR':'H','15N NMR':'N'}
complete = ['HEWL', 'Bovine Beta-Lactoglobulin',
                'Plastocyanin (Anabaena variabilis)',
                'Plastocyanin (Phormidium)',
                'Glutaredoxin', 'CexCD (Apo)',
                'Protein G B1','Xylanase (Bacillus subtilus)']

col=cols[0]
nucl = nuclnames[col]
t = TitrationAnalyser()

#ghost mapping..
DB = PDatabase(server='peat.ucd.ie', username='guest',
               password='123', project='titration_db',
               port=8080)
p=t.extractpKas(DB, col, names=['HEWL'], titratable=False, reliable=False, minspan=0.06)
t.mappKas(DB,col,p,names=['HEWL'],
          nucleus=nucl,calculatespans=False)       




