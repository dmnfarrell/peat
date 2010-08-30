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

shutil.rmtree('peatdb', ignore_errors=True)
path=os.path.abspath('../../..')
peatpath=os.path.abspath('../../../PEATDB')
version = '2.0'

f = Freezer('peatdb', includes=("numpy",),excludes=("wx",))
f.addScript(os.path.join(peatpath, "PEATApp.py"))
f.addScript(os.path.join(peatpath, "Ekin/Ekin_main.py"))
f.addScript(os.path.join(peatpath, "DNAtool/DNAtool.py"))
m=f.mf
f()    # runs the freezing process

'''post freeze'''
#mpl data
import matplotlib
mpldir = matplotlib.get_data_path()
datadir = 'peatdb/mpl-data'
shutil.copytree(mpldir, datadir)

#add peat resource files
resources = ['PEATDB/DNAtool/restriction_enzymes.DAT',            
             'PEATDB/data/AA_masses.txt',
             'PEATDB/App.ico',
             'PEATDB/DNAtool/DNAtool.ico',
             'Protool/AA.DAT',
             'Protool/bbdep02.May.sortlib']
for r in resources:
    shutil.copy(os.path.join(path, r), 'peatdb')
#set icon?

#make zip archive
import zipfile
f = zipfile.ZipFile("peatdb-2.0.zip", "w")
for dirpath, dirnames, filenames in os.walk('peatdb'):
    for fname in filenames:
        fullname = os.path.join(dirpath, fname)        
        f.write(fullname)
f.close()

