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

import string
from Tkinter import *
dict={}
#
#opens a file reads it and prints assignment plus column 2 
#

def find_cols(filename):
    """Get no. of assignment cols in peak file"""
    fd=open(filename)
    lines=fd.readlines()
    fd.close()
    cols=None
    import string
    for line in lines:
        line=string.strip(line)
        if len(line)>0 and 'Assignment' in line:     
            split=line.split()
            print split
            cols = len(split)-1
           
    return cols
    
    
def read_peaks(filename):
    """Get all the peaks"""
    fd=open(filename)
    lines=fd.readlines()
    fd.close()
    peaks={}
    #aas=['D','E','R','K','H',]
    import string
    print '================='
    print 'Reading: %s' %filename
    for line in lines:
        line=string.strip(line)
        if len(line)>6 and line.find('Assignment')==-1:
            if line[0]:
                split=line.split()
                peaks[split[0]]={'N':float(split[1]),'H':float(split[2])}
                #print 'peaks',peaks
                #print 'peaks[split[0]]', peaks[split[0]]                
            else:
                print 'Not reading line',line
    return peaks
            
def calc_delta_shift(H, N, Href, Nref, f=5):  
    """Find combined chem shift"""
    import math
    delta_shift = math.sqrt( math.pow(f*(H - Href),2) + math.pow(N - Nref, 2) )
    return delta_shift
    
def read_all(parentwin,dir,pattern):
    """Read all Sparky peak files in dir that match pattern"""
    import os, re
    files=os.listdir(dir)
    okfiles=[]
    files.sort()
    for file in files:
        if file.find(pattern)!=-1 and file[0]!='.':
            okfiles.append(os.path.join(dir,file))
    print 'files to be used:',okfiles 
    columns=None
    
    # Loop over all files and get peaks    
    dict={}
    regex = "\d+\.\d+"

    for file in okfiles:        
        # Get the pH from the filename        
        filename=os.path.split(file)[-1]
        match = re.findall(regex,filename)
        try:
            pH = float(match[0])
        except:
            'Could not extract the pH from filename'
        
        #Read the file and get the peaks        
        dict[pH]=read_peaks(file)        
    return dict
