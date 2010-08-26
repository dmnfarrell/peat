#
# $Id: pKaSystem_mac_setup.py 211 2007-03-19 12:54:35Z nielsen $
#
# Setup script for pKaTool
#
"""Run: python pKaSystem_mac_setup.py py2app -p Pmw,Pmw.Pmw_1_2,Pmw.Pmw_1_2.lib,Pmw.Pmw_1_2.contrib"""

from distutils.core import setup
import sys
sys.path.append('/usr/local/bin/')
sys.path.append('/Library/Python/2.3/py2app')
import py2app
setup(name='pKaSystem',
      version='1.0',
      description='pKaTool: a package for studying small systems of titratable groups',author='Jens Erik Nielsen',
      author_email='pKa@ucd.ie',
      url='http://enzyme.ucd.ie/Science/pKa/pKaTool',
      license='GNU General Public License v3',
      app=['pKaTool/pKa_system.py'],
      packages=['Pmw','Pmw.Pmw_1_2','Pmw.Pmw_1_2.lib','Pmw.Pmw_1_2.contrib'])
