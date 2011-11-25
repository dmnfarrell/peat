#!/usr/bin/env python
#
# DataPipeline - A data import and fitting tool
# Copyright (C) 2011 Damien Farrell
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact information:
# Email: damien.farrell_at_ucd.ie
# Normal mail:
# Damien Farrell
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland
# 

import os, sys, math, random, numpy, string, types
from datetime import datetime
import ConfigParser, csv
from Base import BaseImporter

"""Custom Importers should be added here."""

class KineticsDataImporter(BaseImporter):
    """This is a custom importer to handle the kinetics data supplied as part of the 
       case study. This data are kinetic assays measured over time intervals per substrate 
       concentration. The data are grouped in rows per time point, each row is represents a
       specific concentration. Each column is one variant.       
       The importer returns a .."""
       
    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):      
        data = {}
        


