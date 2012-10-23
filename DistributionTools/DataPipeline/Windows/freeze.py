#!/usr/bin/env python

#bbfreeze setup file for DataPipeline distribution on Windows
#Damien Farrell, #Nov 2011

"""
This script can be used to create a standalone executable for
either windows or linux. It must be run on the target platform.
You will need to install bbfreeze, see http://pypi.python.org/pypi/bbfreeze/
"""

from bbfreeze import Freezer
import sys, os, shutil

shutil.rmtree('datapipeline', ignore_errors=True)
path=os.path.abspath('../../..')
pipepath=os.path.abspath('../../../DataPipeline')
peatpath=os.path.abspath('../../../PEATDB')
version = '1.2'

f = Freezer('datapipeline', excludes=('wx'))
f.addScript(os.path.join(pipepath, "PipelineApp.py"))
f.addScript(os.path.join(pipepath, "PipelineCommand.py"))
f.addScript(os.path.join(peatpath, "Ekin/ModelDesign.py"))
f.addScript(os.path.join(peatpath, "Ekin/Ekin_main.py"))

#these lines allow the plugins to work
f.addModule('PEATDB.PEATApp')

m=f.mf
f()    # runs the freezing process

'''post freeze'''
#mpl data
import matplotlib
mpldir = matplotlib.get_data_path()
datadir = 'datapipeline/mpl-data'
shutil.copytree(mpldir, datadir)

#add resource files
resources = ['DataPipeline/app.ico',
             'DataPipeline/modeldesign.ico',
             'PEATDB/Ekin/Ekin.ico',
             'PEATDB/Ekin/models.dict']
for r in resources:
    shutil.copy(os.path.join(path, r), 'datapipeline')
tst = 'DataPipeline/testfiles'
shutil.copytree(os.path.join(path, tst), 'datapipeline/testfiles')

#make zip archive
import zipfile
f = zipfile.ZipFile("datapipeline-1.0.zip", "w")
for dirpath, dirnames, filenames in os.walk('datapipeline'):
    for fname in filenames:
        fullname = os.path.join(dirpath, fname)        
        f.write(fullname)
f.close()

