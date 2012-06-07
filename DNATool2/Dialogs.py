#!/usr/bin/env python
#
# DNATool - A program for DNA sequence manipulation
# Copyright (C) 2012- Damien Farrell & Jens Erik Nielsen
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
# Email: farrell.damien_at_gmail.com


"""Dialogs for DNATool"""

from Tkinter import *
import tkSimpleDialog, tkFileDialog, tkMessageBox
import Pmw
import types, os
import tkColorChooser
#import ProgressBar

class TopLevelModalDialog(Toplevel):
    def __init__(self, parent, width=300, height=100):
        Toplevel.__init__(self, parent)
        self.transient(parent)
        self.title('')
        self.parent = parent
        x=(parent.winfo_rootx()+parent.winfo_width()/2)-width/2
        y=(parent.winfo_rooty()+parent.winfo_height()/2)-height/2
        self.geometry('%sx%s+%s+%s' %(width,height,x,y))
        self.body = Frame(self)
        self.initial_focus = self.body
        self.body.pack(fill=BOTH,expand=1,padx=5, pady=5)
        self.grab_set()
        if not self.initial_focus:
            self.initial_focus = self
        self.protocol("WM_DELETE_WINDOW", self.close)
        return

    def close(self):
        self.parent.focus_set()
        self.destroy()
        return

'''class ProgressDialog(TopLevelModalDialog):
    def __init__(self, parent, message='Working', cancel=None):
        TopLevelModalDialog.__init__(self, parent)
        self.parent=parent
        self.title(message)
        progrlbl = Label(self.body,text='Progress:')
        progrlbl.pack(fill=BOTH,padx=2,pady=4)
        self.bar = ProgressBar.ProgressBar(self.body)
        self.bar.frame.pack(fill=Y,padx=2,pady=4)
        if cancel != None:
            self.cancel = Button(self.body,text='cancel',command=cancel)
            self.cancel.pack()

        return

    def updateValue(self, value=None):
        #print 'update'
        #self.after(100, self.updateValue)
        self.update()
        if value!=None:
            self.bar.update(value)
        return'''

class GenericDialog(Frame):
    """Generic dialog"""
    def __init__(self, parent, options=None, defaults=None, callback=None):
        self.parent=parent      
        Frame.__init__(self, parent)
        self.main=self
        if options == None:
            self.options = []
        else:
            self.options = options
        self.defaults = defaults
        self.callback = callback
        self.makeGUI()
        return

    def makeGUI(self):
        """Create GUI from options"""
        self.vars={}
        for item in self.options:
            #print item
            name, widget, values, label = item
            if self.defaults.has_key(name):
                value = self.defaults[name]
            if widget == 'checkbutton':
                var = self.vars[name] = IntVar()
                var.set(value)
                w=Checkbutton(self.main,text=label,variable=var)                
            elif widget == 'entry':
                var = self.vars[name] = StringVar()
                w = Frame(self.main)
                Label(w,text=label+':').pack(side=LEFT,fill=X,expand=1)
                Entry(w, textvariable=var).pack(side=LEFT,fill=X,expand=1)
                var.set(value)
            elif widget == 'scale':
                var = self.vars[name] = IntVar()               
                a,b=values
                w = Frame(self.main)
                Label(w,text=label+':').pack(side=LEFT,fill=X,expand=1)
                Scale(w, from_=a, to=b, orient=HORIZONTAL,
                     variable=var).pack(side=LEFT,expand=1)
                var.set(value)
            elif widget == 'menu':
                var = self.vars[name] = StringVar()
                w = Frame(self.main)
                Label(w,text=label+':').pack(side=LEFT,fill=X,expand=1)
                OptionMenu(w, var, *tuple(values)).pack(side=LEFT,fill=X,expand=1)
                var.set(value)
            elif widget == 'color':
                var = self.vars[name] = StringVar()
                var.set(value)
                def getColor(i):  
                    ctuple, variable = tkColorChooser.askcolor(parent=self, 
                                                            title='Pick a color',
                                                            initialcolor=value)
                    self.vars[i].set(variable)
                w = Button(self.main, text=label, command=lambda: getColor(name))
            w.pack(side=TOP,fill=X,expand=1,padx=2)
        
        self.buttonframe = fr = Frame(self.main)
        fr.pack(side=BOTTOM,fill=BOTH,expand=1)
        b = Button(fr, text="Apply",command=self.callback)
        b.pack(side=LEFT,fill=X,expand=1,padx=1,pady=1)
        return

    def getValue(self, name):
        """Get an option value"""
        if self.vars.has_key(name):
            return self.vars[name].get()

    def getValues(self):
        """Get the current option values as a list of tuples"""
        vals = {}      
        for name in self.vars:
            print name, self.vars[name].get()         
            vals[name] = self.vars[name].get()     
        return vals

    def close(self):
        self.main.destroy()
        self.parent.destroy()
        return

def openFilename(parent, ext=['txt'], savedir=None):
    if not type(ext) is types.ListType:
        ext=[ext]
    if savedir == None:
        savedir = os.getcwd()
    filetypes=[]
    for e in ext:
        filetypes.append(('%s files' %e,'*.%s' %e))
    filetypes.append(("All files","*.*"))
    filename=tkFileDialog.askopenfilename(defaultextension=ext,
                                          initialdir=savedir,
                                          filetypes=filetypes,
                                          parent=parent)
    return filename

def openFilenames(parent, ext=['txt'], savedir=None):
    if not type(ext) is types.ListType:
        ext=[ext]
    if savedir == None:
        savedir = os.getcwd()
    filetypes=[]
    for e in ext:
        filetypes.append(('%s files' %e,'*.%s' %e))
    filetypes.append(("All files","*.*"))
    filename=tkFileDialog.askopenfilenames(defaultextension=ext,
                                          initialdir=savedir,
                                          filetypes=filetypes,
                                          parent=parent)
    return filename

def saveFilename(parent, ext=[''], savedir=None):
    if ext!='':
        filetypes = [('%s files' %ext,'*.%s' %ext)]
    else:
        filetypes = []
    filetypes.append(("All files","*.*"))
    filename=tkFileDialog.asksaveasfilename(defaultextension='.'+ext,
                                          initialdir=savedir,
                                          filetypes=filetypes,
                                          parent=parent)
    return filename

def openDirectory(parent):
    folder = tkFileDialog.askdirectory(parent=parent,
                                        initialdir=os.getcwd(),
                                        title='Select folder')
    return folder



