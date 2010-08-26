#!/usr/bin/env python
#
# pKaTool - analysis of systems of titratable groups
# Copyright (C) 2010 Jens Erik Nielsen
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


import pKa_base
from Tkinter import *
import numpy


class FTIR_data:
    
    def __init__(self,parent):
        
        self.ftir_win = Toplevel()
        self.parent = parent
        self.fit_lines = []
        self.data_points = []
        self.RI_axis_objects = []

        #
        # graph parameters
        #
        self.RIStart=0.0
        self.RIEnd=0.1

        self.pHstart=0.0
        self.pHend=14.0

        self.linewidth=2
        self.data_point_size=5
        
        self.width=900
        self.height=450

        self.graph_startx=80
        self.graph_endx=850
        self.graph_starty=360
        self.graph_endy=60
        #
        # model parameters
        #
        self.offset = DoubleVar()
        self.offset.set(0.0)
        self.c = DoubleVar()
        self.c.set(-1.000)
        # model: ftir = c*sumtits+offset
        
        self.ftir_data = {}
        #
        # draw graph
        #
        self.make_win()

        return

    def fit(self):
        #
        # Set damper
        #
        self.LM_damper = DoubleVar()
        self.LM_damper.set(0.001)
        #
        # Open the window
        #
        self.count_win=Toplevel()
        self.count_win.geometry('100x100+200+200')
        self.count_win.title('FTIR fit in progress')
        self.count=IntVar()
        self.count.set(0)
        Label(self.count_win,text='Iteration #').grid(row=0,column=0) 
        Label(self.count_win,textvariable=self.count).grid(row=0,column=1)
        #
        # Keep track of number of pKa calcs
        #
        self.pka_count=IntVar()
        self.pka_count.set(0)
        Label(self.count_win,text='pKa calc #').grid(row=1,column=0) 
        Label(self.count_win,textvariable=self.pka_count).grid(row=1,column=1)
        #
        # stop fit button
        #
        Button(self.count_win,text='Stop fit',command=self.stop_fit).grid(row=2,column=0,columnspan=2)
        self.count_win.update_idletasks()
        self.keep_running=1
        #
        # Make copies of all the variables
        #
        self.vars=[]
        for group in self.parent.groups:
            self.vars.append(self.parent.groups[group].intpka)
            for group2 in self.parent.groups[group].intenes.keys():
                self.vars.append(self.parent.groups[group].intenes[group2])
        self.vars.append(self.offset)
        self.vars.append(self.c)

        #
        # Start iterating
        #
        print 'Step    Diff   LM_damper'
        for x in range(1,100):
            if not self.keep_running:
                break
            self.fit_LM()
            #
            # Update more stuff
            #
            self.count.set(self.count.get()+1)
            self.count_win.update()
            self.parent.update_scales_from_fit()
            self.parent.titwin.update()
            self.draw_fit()
            #
            # Check convergence
            #
            diff=0
            sum_tits = self.get_sum_of_tits()
            for pH in self.ftir_data.keys():
                diff+=abs(sum_tits[pH]-self.ftir_data[pH])
            
            print '%5d  %6.4f  %6.4f  %6.4f  %6.4f' %(x, diff, self.LM_damper.get(),self.c.get(),self.offset.get())

           # temp = 'Vars: '
           # for v in self.vars:
            #    temp += '%7f  ' %v.get()
            #print temp


            if abs(diff)<0.000001:
                print 'Converged',diff
                break

        self.count_win.destroy()
        self.parent.update_pkasystem_curves()
        return

  
    
    def fit_LM(self):
        """Do Levenberg-Marquardt fitting"""
       
        J,E =self.get_jacobian()
        JT = numpy.transpose(J)
        JTE = numpy.dot(JT,E)
        JTJ = numpy.dot(JT,J)
        JTJd = JTJ + self.LM_damper.get()*numpy.identity(numpy.shape(JTJ)[0])
        invJTJd = numpy.linalg.inv(JTJd)
        q = -numpy.dot(JTE,invJTJd)

        #print 'q',q
        #print 'J',J
        #print 'E',E
        for var in range(len(self.vars)):
            self.vars[var].set(self.vars[var].get()+q[var])
        return

    def stop_fit(self,event=None):
        """Stop the fit"""
        self.keep_running=None
        return

 

    def get_jacobian(self):
        """Get the Jacobian matrix and errors of the data points"""
        #
        # Get the number of data points
        #
        no_data_points = len(self.ftir_data.keys())
        #
        #
        #
        errors = numpy.resize(numpy.array(0,float),[no_data_points])        
        jacobian = numpy.resize(numpy.array(0,float),[no_data_points,len(self.vars)])
        #
        # Precalculate the variation of all parameters
        #
 #       now=self.parent.eval(silent=1)
        now = self.get_sum_of_tits()
        variations=[]

        step = 1e-8

        for var in range(len(self.vars)):
            self.vars[var].set(self.vars[var].get()+step)            
            variations.append(self.get_sum_of_tits(curves=self.parent.eval(silent=1).curves))
            self.vars[var].set(self.vars[var].get()-step)
        #
        # construct jacobian
        #
        data_id=0
     
        sum_of_tits = self.get_sum_of_tits()
        
        for ph in self.ftir_data.keys():                
            errors[data_id] = self.ftir_data[ph]-sum_of_tits[ph]
                
            #
            # Find the derivative at this ph (x) for this data point
            #
            diff=numpy.resize(numpy.array(0,float),[len(self.vars)])
            count=0
            for variation in variations:
                diff[count]=(now[ph]-variation[ph])/step
                count=count+1
            jacobian[data_id]=diff
            data_id=data_id+1
        
        return jacobian,errors
        


    def draw_data(self):
        #
        # draw FTIR data
        #
        for d in self.data_points:
            self.graph.delete(d)
        self.data_points = []
        
        for pH in self.ftir_data.keys():
            x,y = self.get_coordinates(pH, self.ftir_data[pH])
            self.data_points.append(self.graph.create_oval(x,y,x+self.data_point_size,y+self.data_point_size))
        
        return
    
    def draw_fit(self):
        #
        # delete old lines
        #
        for l in self.fit_lines:
            self.graph.delete(l)
        self.fit_lines = []
        #
        # Find min and max values of model and exp. data
        #

        model = self.get_sum_of_tits()
        mina = min(model.values())
        maxa = max(model.values())

        try:
            mine = min(self.ftir_data.values())
            maxe = max(self.ftir_data.values())

            mina = min(mine,mina)
            maxa = max(maxe,maxa)
        except:
            pass
    
        diff = maxa-mina
        maxa+=diff/10.0
        mina-=diff/10.0

        self.RIStart=mina
        self.RIEnd=maxa

        self.draw_RI_axis()
        #
        # draw FTIR fit
        #
        
        phs = model.keys()
        phs.sort()
        for i in range(1, len(phs)):
            x1,y1 = self.get_coordinates(float(phs[i-1]), model[phs[i-1]])
            x2,y2 = self.get_coordinates(float(phs[i]), model[phs[i]])            
            self.fit_lines.append(self.graph.create_line(x1,y1,x2,y2,fill='red',width=self.linewidth))
        

        self.draw_data()
        
        return


    def get_sum_of_tits(self, curves=None):
        #
        # get sum titration curves
        #
        if curves == None:
            curves = self.parent.titration_curves
        
        sum_tits = {}
        for group in curves.keys():
            id=self.parent.ids[group]
            if self.parent.groups[id].active.get()==1:
                for pH in curves[group].keys():
                    if not pH == 'pKa':
                        if sum_tits.has_key(pH):
                            sum_tits[pH]+=curves[group][pH]
                        else:
                            sum_tits[pH]=curves[group][pH]
  
        #write out simulated FTIR curve
  #      phs = sum_tits.keys()
   #     phs.sort()
    #    for i in range(len(phs)):
     #       print '%8f  %8f' %(phs[i], sum_tits[phs[i]]/100)
        
        #apply model
        

        for ph in sum_tits.keys():
            sum_tits[ph] = self.c.get()*sum_tits[ph]+self.offset.get()
        
        
        return sum_tits



    def get_coordinates(self, pH, relative_intensity):
        """ calculates Canvas coordinates for a data point """
        x = self.graph_startx + (pH - self.pHstart)/(self.pHend-self.pHstart)*(self.graph_endx - self.graph_startx)
        y = self.graph_starty - (relative_intensity - self.RIStart)/(self.RIEnd-self.RIStart)*(self.graph_starty-self.graph_endy)
        return x,y
    

    def make_win(self):
        #
        # Buttons
        #
        Label(self.ftir_win,text='Offset').grid(row=0,column=0)
        Entry(self.ftir_win,textvariable=self.offset,width=10).grid(row=0,column=1)
        Label(self.ftir_win,text='Factor').grid(row=0,column=2)
        Entry(self.ftir_win,textvariable=self.c,width=10).grid(row=0,column=3)
        Button(self.ftir_win,text='Update',command=self.draw_fit).grid(row=0,column=4)
        Button(self.ftir_win,text='Reset',command=self.reset).grid(row=0,column=5)
        #
        # Graph
        #
        self.ftir_win.geometry('%dx%d+10+20' %(self.width,self.height))
        self.ftir_win.title('FTIR data')
        self.graph=Canvas(self.ftir_win,bd=5,bg='white',width=self.width,height=self.height,scrollregion=(0,0,self.width,self.height))
        self.graph.grid(row=1,column=0,columnspan=100)

        self.draw_pH_axis()
        self.draw_RI_axis()
        return

    def reset(self):
        self.offset.set(0.0)
        self.c.set(-1.000)
        self.draw_fit()
        return

    def draw_pH_axis(self):
        #
        # pH axis
        #
        self.graph.create_line(self.graph_startx,self.graph_starty,self.graph_endx,self.graph_starty,fill='black',width=self.linewidth)
        self.graph.create_text(self.graph_endx+10,self.graph_starty,text='pH')
        #
        # Tick marks and tick labels
        #
        for tickph in range(int(self.pHstart)*10,int(self.pHend*10),10):
            x,dummy=self.get_coordinates(tickph/10.0,0.0)
            self.graph.create_line(x,self.graph_starty,x,self.graph_starty+5,fill='black',width=self.linewidth)
            self.graph.create_text(x,self.graph_starty+15,text='%5.1f' %(float(tickph)/10.0),fill='black')
        return



    def draw_RI_axis(self):
        
        for o in self.RI_axis_objects:
            self.graph.delete(o)
        self.RI_axis_objects = []
            
        startx,starty=self.get_coordinates(self.pHstart,self.RIStart)
        endx,endy=self.get_coordinates(self.pHend,self.RIEnd)
        #
        # Relative intensity axis
        #
        self.RI_axis_objects.append(self.graph.create_line(startx,starty,startx,endy,fill='black',width=self.linewidth))
        self.RI_axis_objects.append(self.graph.create_text(startx,endy-15,text='Relative intensity / %',fill='black'))
        #
        # Tick marks and tick labels
        #
        start = int(self.RIStart*100000)
        end = int(self.RIEnd*100000)+1
        step = int((end-start)/10.0)
        for tick in range(start,end,step):
            dummy,y=self.get_coordinates(0.0,float(tick)/100000.0)
            self.RI_axis_objects.append(self.graph.create_line(startx,y,startx-5,y,fill='black',width=self.linewidth))
            self.RI_axis_objects.append(self.graph.create_text(startx-25,y,text='%5.4f' %(float(tick)/100000.0),fill='black'))

        return



    def load_ftir_data(self):
        #
        # Load experimental FTIR data
        #
        import tkFileDialog, os
        savedir=os.getcwd()
        filename=tkFileDialog.askopenfilename(defaultextension="*.txt",
                        initialdir=savedir,
                        filetypes=[("Text file","*.txt"),("pKaTool titration curve file","*.tcrv"),("All files","*.*")])

        if not filename:
            return
        self.ftir_data = {}
        
        lines = open(filename).readlines()
        for line in lines:
            if not line[0:1] == '#':
                [ph,value]=line.split()
                self.ftir_data[float(ph)]=float(value)*100
        
        self.parent.update_pkasystem_curves(None,1,None)
        self.draw_fit()
        self.draw_data()
        
        return
                    


