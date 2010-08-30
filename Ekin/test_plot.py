#!/usr/bin/env python

# Written by Damien Farrell, Sep 2007
# Show dataset plot test

import math
import sys,os
from Tkinter import *
from Project import *
import Plotter

sys.path.append('/software/python/')

class testplotter(Frame,Plotter.PEAT_Plot):
    """Test plotter for doing ekin plots"""
    def __init__(self,parent=None):
        
        Frame.__init__(self,parent=None,width=400, height=400)

        import tkFont
        self.detailsfont = tkFont.Font ( family="Helvetica", size=12, weight="bold" )
    
        #for storing all currently loaded proteins and their datasets    
        self.proteindatasets=[]    
        row=0
        column=0
        #load the data for the plots and dataset selector
        self.get_Ekin_data()
        #
        #Do widget stuff now
        #label
        Label(self,text='test',font=self.detailsfont).grid(row=row,column=column,sticky='NWS')
         
        #
        # Create the canvas and widgets needed to plot
        #            
        self.gcv=Canvas(self,bd=5,bg='white',
                           width=420,
                           height=420)
        self.gcv.grid(row=row+1,column=0,sticky='news',padx=6,pady=4)               
        self.rowconfigure(1,weight=1)
        self.columnconfigure(0,weight=1)
        #callback for resizing
        self.gcv.bind("<Configure>", self.checksize)
        
        #
        # Add the menu for choosing the datasets, this bit taken from Ekin_main
        #
        self.selected_dataset=StringVar()
        self.selected_dataset.set('None')
        self.dataset_select_button=Menubutton(self,textvariable=self.selected_dataset,relief=RAISED)
        self.dataset_sel_menu=Menu(self.dataset_select_button,tearoff=0)
        self.dataset_select_button['menu']=self.dataset_sel_menu
            
        self.ekinplotter=Plotter.PEAT_Plot(self,canvas=self.gcv)
        self.ekinplotter.plot_ekin_data(self.plotdata,self.fit_data)
        
        return        

    def checksize(self,event):
        """Call when widget resized"""
        #global w,h,xp,yp,xdelta,ydelta
        w = event.width
        h = event.height
        self.ekinplotter.rescale_canvas(w,h)
        return
         
    def get_Ekin_data(self):
        """Open the Ekin data from an ekin proj file"""
        ekinfilename='sampled66a.Ekinprj'        
        ekinfile=os.getcwd()
        print 'opened Ekin project: '+ekinfilename
        P = Project()
        self.plotdata,self.fit_data = P.open_project_nogui(filename=os.path.join(ekinfile,ekinfilename))        
        #print self.plotdata
        #print self.fit_data
        return
        
class Container(Frame):
     def __init__(self, parent=None):
         Frame.__init__(self, parent)
         self.pack()
         X = testplotter(self)
         X.pack(fill=BOTH, expand=YES)
         
         
if __name__=="__main__":
    Container().mainloop()
    

