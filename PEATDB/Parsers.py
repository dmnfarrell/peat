#!/usr/bin/env python
#
# Protein Engineering Analysis Tool DataBase (PEATDB)
# Copyright (C) 2010 Damien Farrell & Jens Erik Nielsen
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
# Email: Jens.Nielsen_at_gmail.com 
# Normal mail:
# Jens Nielsen
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland
# 
# Author: Damien Farrell 2011

import urllib
import sys
import xml.parsers.expat as expat

class XMLParser(object):
    def __init__(self):
        self._parser = expat.ParserCreate()
        self._parser.StartElementHandler = self.start
        self._parser.EndElementHandler = self.end
        self._parser.CharacterDataHandler = self.data

    def feed(self, data):
        self._parser.Parse(data, 0)
        

    def close(self):
        self._parser.Parse("", 1) # end of data
        del self._parser # get rid of circular references

    def start(self, tag, attrs):
        #print "START", repr(tag), attrs
        self.current = {tag: attrs}
   
    def end(self, tag):
        #print "END", repr(tag)
        return

    def data(self, data):
        #print "DATA", repr(data)
        return

    def openurl(self, url):
        u = urllib.urlopen(url)    
        data = u.read()
        self.feed(data)
        self.close()
        return self.current
    
if __name__ == '__main__':    
    url = 'http://www.rcsb.org/pdb/rest/describePDB?structureId=4hhb'
    p = XMLParser()
    d = p.openurl(url)
    for k in d['PDB']:
        print k, d['PDB'][k]
