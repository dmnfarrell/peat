#!/usr/bin/env python

#bbfreeze setup file for PEAT_DB distribution on Windows
#Damien Farrell, #October 2009

"""
This script can be used to create a standalone executable for
either windows or linux. It must be run on the target platform.
You will need to install bbfreeze, see http://pypi.python.org/pypi/bbfreeze/
"""

from bbfreeze import Freezer
import sys, os, shutil

shutil.rmtree('DNAtool', ignore_errors=True)
path=os.path.abspath('../../..')
dnatoolpath=os.path.abspath('../../../PEATDB/DNAtool')
version = '2.0'

f = Freezer('DNAtool', excludes=('wx'))
f.addScript(os.path.join(dnatoolpath, "DNAtool.py"))
#these lines allow the plugins to work

m=f.mf
f()    # runs the freezing process

'''post freeze'''

#add resource files
resources = ['restriction_enzymes.DAT',
             'DNAtool.ico',
             'test.DTP','DNAseq.seq']

for r in resources:
    shutil.copy(os.path.join(dnatoolpath, r), 'DNAtool')

#make zip archive
import zipfile
f = zipfile.ZipFile("DNAtool-2.0.zip", "w")
for dirpath, dirnames, filenames in os.walk('peatdb'):
    for fname in filenames:
        fullname = os.path.join(dirpath, fname)        
        f.write(fullname)
f.close()

