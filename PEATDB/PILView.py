#!/usr/bin/env python
#
# This file is part Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-
# All rights reserved
#
#
# Written by D Farrell, Sep 2008
#

from Tkinter import *
import tkFileDialog, tkMessageBox, tkSimpleDialog
import Pmw 
import re
import os
import time
import math
try:
    from PIL import Image, ImageTk
except:
    print 'You need PIL installed. If using fedora do yum install python-imaging*'
    print
from Prefs import Preferences

class PILViewer(Frame):
    """
    Basic image viewing app using PIL for python
    """    
    appname = 'PIL viewer'

    def __init__(self, parent=None, peatinfo=None, title=None,
                     data=None, imgfile=None):
        "Initialize the application."
        self.parent=parent
        #If there is data to be loaded, show the dialog first
        if not self.parent:
            Frame.__init__(self)
            self.pilwin=self.master
            self.peatinfo=None
        else:
            self.pilwin=Toplevel()
            #self.peatinfo=peatinfo      #reference to peat protein/field
        if title != None:
            self.title = 'imgviewer_' + title
        else:
            self.title = 'PIL Image Viewer'
        self.pilwin.title(self.title)
        import platform
        self.currplatform=platform.system()        
        if not hasattr(self,'defaultsavedir'):
            self.defaultsavedir = os.getcwd()
        self.preferences=Preferences('PILViewer',{'check_for_update':1})

        self.pilwin.geometry('+200+100')
        self.currentimage = None
        self.setupGUI()
        
        if imgfile != None:
            self.openFile(imgfile)
        elif data != None:
            self.currentimage =  Image.frombuffer('RGB', (100,100), data, "raw", 'RGB', 0, 1) 
            self.updateCanvas()
        return

    def setupGUI(self):
        """Do GUI elements"""
        self.createMenuBar()
        self.createCanvas()
        self.pilwin.bind("<Configure>", self.resize)        
        self.pilwin.rowconfigure(0,weight=1)
        self.pilwin.columnconfigure(0,weight=1)
        return
        
    def createMenuBar(self):
        """Create the menu bar for the application. """
        self.menu=Menu(self.pilwin)
        self.file_menu={                      
                        '01Quit':{'cmd':self.quit}}
        
        self.file_menu = {'01Open Image':{'cmd':self.openFile},
                          '02Save As':{'cmd':self.saveFile}}
        self.file_menu=self.create_pulldown(self.menu,self.file_menu)
        self.menu.add_cascade(label='File',menu=self.file_menu['var'])

                     
        self.help_menu={'01Online Help':{'cmd':self.online_documentation},
                        '02About':{'cmd':self.about_PIL}}        
        self.help_menu=self.create_pulldown(self.menu,self.help_menu)
        self.menu.add_cascade(label='Help',menu=self.help_menu['var'])
        
        self.pilwin.config(menu=self.menu)
     
        return
    
    def create_pulldown(self,menu,dict):
        #
        # Create a pulldown in var from the info in dict
        #
        var=Menu(menu,tearoff=0)
        items=dict.keys()
        items.sort()
        for item in items:
            if item[-3:]=='sep':
                var.add_separator()
            else:
                #
                # Do we have a command?
                #
                command=None
                if dict[item].has_key('cmd'):
                    command=dict[item]['cmd']
                #
                # Put the command in there
                #
                if dict[item].has_key('sc'):
                    var.add_command(label='%-25s %9s' %(item[2:],dict[item]['sc']),command=command)
                else:
                    var.add_command(label='%-25s' %(item[2:]),command=command)
        dict['var']=var
        return dict

    def createCanvas(self):
        """Create canvas for images"""
        self.imgcanvas = Canvas(self.pilwin, width=600, height=400, bg='gray25')
        self.imgcanvas.grid(row=0,column=0,sticky='news')
        canv = self.imgcanvas
        self.Yscrollbar = AutoScrollbar(self.pilwin,orient=VERTICAL,command=self.set_yviews)
        self.Yscrollbar.grid(row=0,column=1,rowspan=1,sticky='news',pady=0,ipady=0)
        self.Xscrollbar = AutoScrollbar(self.pilwin,orient=HORIZONTAL,command=self.set_xviews)
        self.Xscrollbar.grid(row=1,column=0,columnspan=1,sticky='news')
        canv['xscrollcommand'] = self.Xscrollbar.set
        canv['yscrollcommand'] = self.Yscrollbar.set        
        return

    def set_xviews(self,*args):
        """Set the xview of table and col header"""
        apply(self.imgcanvas.xview,args)        
        return
        
    def set_yviews(self,*args):
        """Set the xview of table and row header"""
        apply(self.imgcanvas.yview,args)        
        return   

    def displayFile(self, img, width=None, height=None):
        """Show the image file on the canvas"""
        self.imgcanvas.delete('img')
        top=0;left=0
        if width != None:
            fw=int(width)
            fh=int(height)            
            img, top, left = self.doImageRescale(img, fw, fh)
        self.tkimg = ImageTk.PhotoImage(img)
        imgobj = self.imgcanvas.create_image(0, 0, image=self.tkimg, anchor='nw', tag='img')
        
        return

    def resize(self, event=None):
        #self.pilwin.after(200)
        #if time
        self.updateCanvas(event)
        
    def updateCanvas(self, event=None):
        """We want the image to fit to the window"""
        if event != None:
            w = event.width
            h = event.height
        else:
            w=None;h=None       
        if self.currentimage != None:
            self.displayFile(self.currentimage, width=w, height=h)
        return    

    def doImageRescale(self, img, framew, frameh):
        """ensures that the image is displayed with the
           correct ratio as the frame is resized """   
        
        imgw, imgh = img.size
        imageratio = float(imgw)/float(imgh)
        print imageratio
        if (imageratio >= 1):
            w = framew
            h = int(w / imageratio)
            left=0
            top = frameh/2 - h/2            
            if (h > frameh):
                h = frameh
                w = int(h * imageratio)
                top=0
                left=framew/2 - w/2 
        elif (imageratio < 1):
            h = frameh
            w = int(h * imageratio)
            if (w > framew):
                w = framew
                h = int(w / imageratio)
            top=0    
            left = framew/2 - w/2 
        newimg = img.resize((w, h), Image.ANTIALIAS)        
        return newimg, top, left        
    
    def openFile(self, filename=None):
        import os
        if filename == None:
            import tkFileDialog 
            filename=tkFileDialog.askopenfilename(defaultextension='.jpg',
                                                      initialdir=os.getcwd(),
                                                      filetypes=[("jpg file","*.jpg"),
                                                                 ("png file","*.png"), 
                                                                 ("All files","*.*")],
                                                      parent=self.pilwin)
        if os.path.isfile(filename):
            self.filename = filename
            self.currentimage = Image.open(filename)
            self.updateCanvas()

        return

    def saveFile(self):
        if self.currentimage != None:
            import tkFileDialog 
            filename=tkFileDialog.asksaveasfilename(defaultextension='.jpg',
                                                      initialdir=os.getcwd(),
                                                      filetypes=[("jpg file","*.jpg"),
                                                                 ("png file","*.png"), 
                                                                 ("All files","*.*")],
                                                      parent=self.pilwin)
            if filename:
                self.currentimage.save(filename)

        return
    
    def online_documentation(self):

        return
    
    def about_PIL(self):

        return


class AutoScrollbar(Scrollbar):
    # a scrollbar that hides itself if it's not needed.  only
    # works if you use the grid geometry manager.
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            # grid_remove is currently missing from Tkinter!
            self.tk.call("grid", "remove", self)
        else:
            self.grid()
        Scrollbar.set(self, lo, hi)
    def pack(self, **kw):
        raise TclError, "cannot use pack with this widget"
    def place(self, **kw):
        raise TclError, "cannot use place with this widget"
    
def main():
    "Run the application."
    
    import sys, os  
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="imgfile", 
                        help="Open a table file", metavar="FILE")    
    opts, remainder = parser.parse_args()
    if opts.imgfile != None:
        app=PILViewer(imgfile=opts.imgfile)
    else:    
        app=PILViewer()
    app.mainloop()
    return

if __name__ == '__main__':
    main()
