from setuptools import setup
import sys,os
home=os.path.expanduser('~')

setup(name='PEATDB', version='2.2',
	description='Protein Engineering Analysis Tool',
	author='Nielsen Group, UCD',
	author_email='farrell.damien[at]gmail.com',
	url='http://enzyme.ucd.ie/PEAT',
	license='GPL v3',
	packages=['PEATDB','PEATDB.DNAtool','PEATDB.Ekin','PEATDB.plugins',
		    'PEATSA','PEATSA.Core','PEATSA.WebApp',
		    'Protool','pKaTool'],
	package_data={'PEATDB.DNAtool':['restriction_enzymes.DAT','DNAtool.ico'], 
		      'PEATDB':['logo.gif','App.ico','data/AA_masses.txt'],
	              'PEATDB.Ekin':['models.dict'],
		      'Protool':['AA.DAT','bbdep02.May.sortlib']},
    install_requires=['numpy>=1.1',
                      'matplotlib>=0.98.5.2','ZODB3',
                      'MySQL-python','relstorage>=1.4'],
    entry_points = { 'gui_scripts': [
                     'peatdb = PEATDB.PEATApp:main',
                     'ekin= PEATDB.Ekin.Ekin_main:main'],}
)
