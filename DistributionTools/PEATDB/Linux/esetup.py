#more advanced version of the distutils setup.py using setuptools

from setuptools import setup, find_packages
import ez_setup
ez_setup.use_setuptools()

setup(
    name='PEAT_DB', version='2.0',
    description='Protein Engineering Analysis Tool',
    author='Nielsen Group, UCD',
    author_email='PEAT_support@ucd.ie',
    url='http://enzyme.ucd.ie/PEAT',
    license='For academic and non-profit use only. No redistribution, no alterations.',
    keywords = 'peat database protein bio',
    packages=['PEAT_DB','PEAT_DB.DNAtool','PEAT_DB.Ekin','PEAT_DB.Search',
              'PEAT_DB.Search.protein_properties','Protool','pKaTool'],
    package_dir={'':'/var/www/python'},
    package_data={'DNAtool':['restriction_enzymes.DAT','DNAtool.ico'],
                  'PEAT_DB':['logo.gif','DB_Main.ico','data/AA_masses.txt'],
                  'Protool':['AA.DAT','bbdep02.May.sortlib']},
    install_requires = ['Pmw>=1.2','numpy>=1.1','matplotlib>=0.98.5'], #'pysvn>=1.6.0'
    #dependency_links = [ "http://peak.telecommunity.com/snapshots/"],
    entry_points = { 'gui_scripts': [ 'peat = PEAT_DB.DB_Main'] }
)
