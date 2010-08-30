#!/bin/env python
#
# This file is part Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-
# All rights reserved
#
# Damien Farrell April 2010

from Tkinter import *
import Pmw

class YasaraControl(Frame):
    """A yasara controller for comparison use"""
    def __init__(self, parent, yasara):        
        Frame.__init__(self, parent, height=200,width=160)
        self.yasara = yasara
        c=Button(self,text='Align', command=self.align)
        c.grid(row=1,column=0,sticky='news',padx=2,pady=2)        
        '''c = Pmw.Counter(parent,
                labelpos = 'w',
                label_text = 'residue:',
                entryfield_value = 0,
                entryfield_command = self.selectRes,
                entryfield_validate = {'validator' : 'integer',
                            'min' : 0, 'max' : 1000})
        c.grid(row=2,column=0,columnspan=2,padx=2,pady=2)'''
       
        return

    def align(self):
        """try to align objects"""
        Y = self.yasara
        Y.AlignMultiAll()
        return

    def selectRes(self):
        """Allow highlight residue from list"""
        
        return
