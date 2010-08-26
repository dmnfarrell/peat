#
# $Id: pKaTool_mac_setup.py 211 2007-03-19 12:54:35Z nielsen $
#
# Setup script for pKaTool
#
"""Run: python pKaTool_mac_setup.py py2app -p Pmw,Pmw.Pmw_1_2,Pmw.Pmw_1_2.lib,Pmw.Pmw_1_2.contrib"""

from distutils.core import setup
import sys
sys.path.append('/usr/local/bin/')
import py2app
setup(name='pKaTool',
      version='1.0',
      description='pKaTool: a package for analysing the output from protein pKa calculations and for studying small systems of titratable groups',author='Jens Erik Nielsen',
      author_email='pKa@ucd.ie',
      url='http://enzyme.ucd.ie/Science/pKa',
      license='GNU General Public License v3',
      app=['pKaTool/pKaTool_main.py'],
      packages=['Pmw','Pmw.Pmw_1_2','Pmw.Pmw_1_2.lib','Pmw.Pmw_1_2.contrib'])
