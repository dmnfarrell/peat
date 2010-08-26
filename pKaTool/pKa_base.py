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

try:
    from Tkinter import *
except:
    pass

import sys,os

class pKa_base:

    def set_graph_size(self,x_max,x_min,y_max,y_min,x_tick_level=1.0,y_tick_level=0.25,width=1200,height=450):
        self.y_max=y_max
        self.y_min=y_min
        self.x_max=x_max
        self.x_min=x_min
        self.graph_y_size=height-150
        self.graph_x_size=width-150
        self.fact_y=self.graph_y_size/(self.y_max-self.y_min)
        self.fact_x=self.graph_x_size/(self.x_max-self.x_min)
        self.x_tick_level=x_tick_level
        self.y_tick_level=y_tick_level
        return

    #
    # ----
    #

    def get_xy(self,pH,crg):
        #
        # Get x and y coordinates
        #
        if not hasattr(self,'y_min'):
            self.set_graph_size(20,0.00001,1.0,-1.0)
        y=crg-self.y_min
        y=50+self.graph_y_size-self.fact_y*y
        return self.get_x(pH),y

    #
    # ----
    #

    def get_x(self,pH):
        #
        # pH
        #
        return (pH-self.x_min)*self.fact_x+70

    #
    # -----------------
    #

    def draw_ordinates(self,canvas,y_label='Charge',x_label='pH'):
        #
        # Draw the axes
        #
        if not hasattr(self,'y_min'):
            self.set_graph_size(20,0.00001,1.0,-1.0)
        startx,starty=self.get_xy(self.x_min,self.y_min)
        endx,endy=self.get_xy(self.x_max,self.y_max)
        y_axis=starty
        x_axis=startx
        #
        # pH axis
        #
        self.pH_axis(canvas,startx,starty,endx,endy,x_label)
        #
        # Charge axis
        #
        canvas.create_line(x_axis,y_axis,x_axis,endy,fill='black',width=self.linewidth)
        canvas.create_text(x_axis,endy-15,text=y_label,fill='black')
        #
        # Tick marks and tick labels
        #
        for tickcrg in range(int(self.y_min*100),int(self.y_max*100)+int(self.y_tick_level*100.0),int(self.y_tick_level*100.0)):
            dummy,y=self.get_xy(0.0,float(tickcrg)/100.0)
            canvas.create_line(x_axis,y,x_axis-5,y,fill='black',width=self.linewidth)
            canvas.create_text(x_axis-25,y,text='%5.2f' %(float(tickcrg)/100.0),fill='black')
        return
    
    #
    # ----------------
    #
    
    def pH_axis(self,canvas,startx,starty,endx,endy,x_label='pH'):
        #
        # Draw the pH axis
        #
        y_axis=starty
        x_axis=startx
        #
        # pH axis
        #
        canvas.create_line(startx,starty,endx,starty,fill='black',width=self.linewidth)
        canvas.create_text(endx+10,starty,text=x_label)
        #
        # Tick marks and tick labels
        #
        for tickph in range(int(self.x_min)*10,int(self.x_max*10.0)+int(self.x_tick_level*10.0),int(self.x_tick_level*10.0)):
            x,dummy=self.get_xy(float(tickph/10.0),self.y_min)
            canvas.create_line(x,y_axis,x,y_axis+5,fill='black',width=self.linewidth)
            canvas.create_text(x,y_axis+15,text='%5.1f' %(float(tickph)/10.0),fill='black')
        return

    #
    # ----------
    #
    
    def update_curves(self,curves=None,labels=None):
        #
        # Delete all lines from last round
        #
        if not getattr(self,'lines',None):
            self.lines={}
        #
        for line in self.lines.keys():
            self.tc.delete(line)
            del self.lines[line]
        #
        # If no curves then return
        #
        if not curves:
            return
        #
        # Define colour_order if it's not defined already
        #
        if not getattr(self,'colour_order',None):
            self.colour_order=['black','red','blue','green','grey','magenta','cyan']
        #
        # Draw the titration curves
        #
        groups=curves.keys()
        groups.sort()
        #pHvalues=curves[groups[0]].keys()
        #pHvalues.sort()
        #pHvalues=pHvalues[:-1]
        #
        # Keep track of the label positions
        #
        self.label_pos=[]
        done_label={}
        #
        # Start looping
        #
        group_count=0
        for group in groups:
            pHvalues=curves[group].keys()
            pHvalues.sort()
            pHvalues=pHvalues[:-1]
            lastpH=pHvalues[0]
            lastcrg=curves[group][lastpH]
            colour=self.colour_order[group_count%len(self.colour_order)]
            #
            #
            pkadone=None
            for pH in pHvalues[1:]:
                lastx,lasty=self.get_xy(lastpH,lastcrg)
                crg=curves[group][pH]
                x,y=self.get_xy(pH,crg)
                self.lines[(self.tc.create_line(lastx,lasty,float(x),float(y),
                                                      fill=colour,
                                                      width=self.linewidth))]=1
                lastcrg=crg
                lastpH=pH
                #
                # Label and pKa value
                #
                pKa=curves[group]['pKa']
                if abs(pKa-pH)<=abs(pHvalues[0]-pHvalues[1]) and not pkadone:
                    ty=float(y)#+float(group_count)*10.0
                    newx,newy=self.resolve_clash(x,ty)
                    #if newx!=x or newy!=ty:
                    #    #
                    #    # Draw a line to the label position
                    #    #
                    #    self.lines[self.tc.create_line(x,y,newx-50,newy,fill=colour,
                    #                                   width=self.linewidth)]=1
                    x=self.tc.create_text(newx,newy,text='%s, pKa: %.1f' %(group,pKa),fill=colour)
                    self.lines[x]=1
                    self.label_pos.append([newx,newy])
                    pkadone=1
            #
            # Did we do a label for this group?
            #
            if not pkadone:
                x,y=self.get_xy(pHvalues[0],curves[group][pHvalues[0]])
                x=x+50
                y=y-10
                newx,newy=self.resolve_clash(x,y)
                if newx!=x or newy!=y:
                    #
                    # Draw a line to the label position
                    #
                    self.lines[self.tc.create_line(x,y,newx,newy,fill=colour,
                                                   width=self.linewidth)]=1
                x=self.tc.create_text(newx,newy,text='%s, pKa: %.1f' %(group,pKa),fill=colour)
                self.lines[x]=1
                self.label_pos.append([newx,newy])
            #
            # Update the counter for colours
            #
            group_count=group_count+1
        return

    #
    # -----
    #

    def resolve_clash(self,x,y):
        """Resolve label clashes if any..."""
        label_size_x=100
        label_size_y=20
        orgy=y
        orgx=x
        clash_x=0
        clash_y=0
        counter=0
        first_check=1
        while (clash_x>0 or clash_y>0 or first_check) and counter<200 :
            clash_x=0
            clash_y=0
            for xold,yold in self.label_pos:
                diffx=abs(x-xold)
                diffy=abs(y-yold)
                if diffx<label_size_x and diffy<label_size_y:
                    # Clash
                    #
                    # Record the smallest clash distance
                    if abs(label_size_x-diffx) < abs(label_size_y-diffy):
                        clash_x=clash_x+1
                    else:
                        clash_y=clash_y+1
            #
            # Resolve clash
            #
            counter=counter+1
            if not first_check:
                if clash_y>0:
                    y=y+10
                    if y>400:
                        x=x+10
                        y=orgy
                if clash_x>0:
                    x=x+10
            else:
                first_check=None
        #
        # return the ok positions
        #
        return x,y
    #
    # ------
    #
    

    def menubutton_list(self,window=None,variable=None,list=None,default=None,indicatoron=0):
        variable.set(default)
        # the button
        button=Menubutton(window,textvariable=variable,relief=RAISED)
        # the menu
        menu=Menu(button,tearoff=0)
        button['menu']=menu
        for type in list:
            menu.add_radiobutton(label=type,variable=variable,value=type,indicatoron=indicatoron)
        return button
