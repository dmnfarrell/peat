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

"""Command line script for data pipeline."""
from Base import Pipeline
import os, random

def loadProject(filename):
    import pickle
    try:
        f = open(filename,'r')
        p = pickle.load(f)
        return p
    except Exception,e:
        print 'failed to load project'
        print 'Error returned:', e
        return

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-c", "--conf", dest="conf",
                            help="Provide a conf file", metavar="FILE")
    parser.add_option("-f", "--file", dest="file",
                        help="Raw file", metavar="FILE")
    parser.add_option("-d", "--dir", dest="directory",
                        help="Folder of raw files")    
    parser.add_option("-p", "--project", dest="project",
                        help="Project file", metavar="FILE")

    opts, remainder = parser.parse_args()
    P = Pipeline()
    if opts.project != None:
        P = loadProject(opts.project)
    else:
        if opts.conf != None:
            P.parseConfig(opts.conf)
        if opts.file != None:
            P.openRaw(opts.file)
        if opts.directory != None:
            P.addFolder(opts.directory)
    P.run()
    
if __name__ == '__main__':
    main()
