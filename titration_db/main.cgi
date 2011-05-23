#!/bin/env python
#
#
# This file is part Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-
# All rights reserved
#
#
# Rewritten by D Farrell, Jan 2009
#

import cgi
import cgitb
cgitb.enable()
import sys,os
sys.path.append('/var/www/python')
#sys.path.append('/home/people/farrell/python')
#from WebInterface import PEATWeb
#from PEATDB.web import PEATWeb
from titdb import titdbWeb

#Execute display of PEAT project
P=titdbWeb(server='peat.ucd.ie',
    	  project='titration_db',
          user='guest',
          passwd='123',
	  port=8080,
          bindir='/cgi-bin/titration_db',
          fullpath='/var/www/html/titration_db')
