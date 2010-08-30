#more advanced version of the distutils setup.py using setuptools

from setuptools import setup, find_packages
#import ez_setup
#ez_setup.use_setuptools()
import sys, os
path=os.path.abspath('../../..')

setup(
    name='PEATDB', version='2.0',
    description='Protein Engineering Analysis Tool',
    author='Nielsen Group, UCD',
    author_email='PEAT_support@ucd.ie',
    url='http://enzyme.ucd.ie/PEAT',
    license='For academic and non-profit use only. No redistribution, no alterations.',
    keywords = 'peat database protein bio',
    packages=['PEATDB','PEATDB.DNAtool','PEATDB.Ekin','Protool'],
    package_dir={'': path},
    package_data={'DNAtool':['restriction_enzymes.DAT','DNAtool.ico'],
                  'PEATDB':['logo.gif','App.ico','data/AA_masses.txt'],
                  'Protool':['AA.DAT','bbdep02.May.sortlib']},
    install_requires = ['numpy','matplotlib','zodb3'],  
    entry_points = { 'gui_scripts': [ 'peat = PEATDB.PEATApp'] }
)
