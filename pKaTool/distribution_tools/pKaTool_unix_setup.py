#
# $Id: pKaTool_setup.py 210 2007-03-18 14:30:51Z nielsen $
#
# Setup script for pKaTool
#
from distutils.core import setup
import os
        
setup(name='pKaTool',
      version='1.1',
      description='pKaTool: a package for analysing protein pKa calculations, small systems of titratable groups and for fitting NMR titration data',author='Jens Erik Nielsen',
      author_email='pKa@ucd.ie',
      url='http://enzyme.ucd.ie/Science/pKa/pKaTool',
      license='GNU General Public License v3',
      packages=['pKaTool','pKarun','Protool','Pmw','Pmw.Pmw_1_2','Pmw.Pmw_1_2.lib','Pmw.Pmw_1_2.contrib'],
      package_data={'pKaTool':['pKaTool.ico'],'Protool':['bbdep02.May.sortlib']})
