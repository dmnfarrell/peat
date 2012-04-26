#more advanced version of the distutils setup.py using setuptools

from setuptools import setup, find_packages
import sys, os
path=os.path.abspath('../../..')

setup(
    name='PEATDB', version='2.2',
    description='Protein Engineering Analysis Tool',
    author='Nielsen Group, UCD',
    author_email='farrell.damien[at]gmail.com',
    url='http://enzyme.ucd.ie/PEAT',
    license='GPL v3',
    keywords = 'peat database protein biology',
    packages=['PEATDB','PEATDB.DNAtool','PEATDB.Ekin','Protool'],
    package_dir={'': path},
    package_data={'DNAtool':['restriction_enzymes.DAT','DNAtool.ico'],
                  'PEATDB':['logo.gif','App.ico','data/AA_masses.txt'],
                  'Protool':['AA.DAT','bbdep02.May.sortlib']},
    install_requires = ['numpy','matplotlib','zodb3','MySQLdb','relstorage'],  
    entry_points = { 'gui_scripts': [ 'peat = PEATDB.PEATApp'] }
)
