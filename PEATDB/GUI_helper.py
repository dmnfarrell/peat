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

from Tkinter import *

class GUI_help:

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
                # Do we have a command?
                command=None
                if dict[item].has_key('cmd'):
                    command=dict[item]['cmd']
                # Put the command in there
                if dict[item].has_key('sc'):
                    var.add_command(label='%-25s %9s' %(item[2:],dict[item]['sc']),command=command)
                else:
                    var.add_command(label='%-25s' %(item[2:]),command=command)

        dict['var']=var
        return dict

    def set_pulldown_and_shortcut_states(self,states,controlvar):
        """
        # Change the states of pulldown menu items
        # states is a dictionary specifying the new states of each menu
        # controlvar is the control variable that specified the menu in the first place
        """
        var=controlvar['var']
        import string
        for item in states.keys():
            found=False
            for menuitem in controlvar.keys():
                #
                #
                text=menuitem[2:]
                if text==item:
                    found=True
                    # Get the menu number
                    number=int(menuitem[:2])-1
                    # Get the keyboard shortcut
                    sc=None
                    if controlvar[menuitem].has_key('sc'):
                        sc=controlvar[menuitem]['sc']
                    if sc:
                        sc_full=self.get_control_seq(sc) # Expand the control sequence
                    #
                    # Set the states and bindings
                    #
                    if states[item]:
                        var.entryconfigure(number,state=NORMAL)
                        if sc:
                            self.master.bind('<%s>' %sc_full,controlvar[menuitem]['cmd']) # Bind
                    else:
                        var.entryconfigure(number,state=DISABLED)
                        if sc:
                            self.master.unbind('<%s>' %sc_full)
                    break
            if not found:
                print 'Could not update state for',item
        self.master.update()
        return

    def get_control_seq(self,shorthand):
        #
        # Expand a shorthand control sequence (e.g. Ctrl+f)
        # to a Tkinter sequence (e.g. Control-KeyPress-f)
        #
        import string
        split=string.split(shorthand,'+')
        dict={'Ctrl':'Control','Del':'Delete'}
        try:
            long=dict[split[0]]
        except KeyError:
            return None
        #
        # Put the key info in there
        #
        if len(split)==2:
            long=long+'-KeyPress-'+string.lower(split[1])
        return long

    def get_text_input(self,parent,fields,title=None):
        """
        # Open a window and get text input from the user
        # The fields to request are in the fields list
        """
        if not getattr(self,'window',None):
            self.window=[]
        if not getattr(self,'winvars',None):
            self.winvars=[]
        #
        # OPen window
        #
        thisid=Toplevel()
        thisid.geometry('+%d+%d' %(parent.winfo_rootx()+50,
                                   parent.winfo_rooty()+50))
        self.window=thisid

        efields=[]
        row=1
        if title:
            Label(thisid,text=title).grid(row=row,column=0)
            row=row+1
        for field,size in fields:
            self.winvars.append(StringVar())
            Label(thisid,text=field).grid(row=row,column=0,sticky='w')
            #
            E=Entry(thisid,textvariable=self.winvars[-1],width=size)
            E.grid(row=row,column=1,sticky='w')
            E.bind("<Return>", self.close_thiswin)
            efields.append(E)
            row=row+1

        Button(thisid,text='Cancel',command=self.cancel_thiswin).grid(row=row,column=0)
        Button(thisid,text='Ok',command=self.close_thiswin).grid(row=row,column=1,sticky='news')
        #
        # Bindings
        #
        thisid.bind("&lt;Return>", self.close_thiswin)
        #
        # Grab the focus
        #
        thisid.lift()
        thisid.grab_set()
        efields[0].focus_set()
        parent.wait_window(thisid)
        results=[]
        if self.winvars:
            for var in self.winvars:
                results.append(var.get())
            self.winvars=[]
            return results
        return None

    def close_thiswin(self,event=None):
        self.window.destroy()
        return

    def cancel_thiswin(self):
        self.winvars=None
        self.window.destroy()
        return

    def set_geometry(self,pwidget,widget):
        """Set the position of a widget in the middle of pwidget"""
        w,h,x,y=self.get_geometry(pwidget)
        sw,sh,dummy,dummy2=self.get_geometry(widget)
        xoffset=int((w-sw)/2)
        yoffset=int((h-sh)/2)
        widget.geometry('+%d+%d' %(x+xoffset,y+yoffset))
        return

    def get_geometry(self, widget):
        """Get the geometry of a widget
        Return width,height,xorg,yorg"""
        widget.update_idletasks()
        txt=widget.winfo_geometry()
        width=int(txt.split('x')[0])
        rest=txt.split('x')[1]
        height=int(rest.split('+')[0])
        xorg=int(rest.split('+')[1])
        yorg=int(rest.split('+')[2])
        return width,height,xorg,yorg

    def set_centered_geometry(self,parent,widget):
        """Place a widget on top of a parent"""
        w,h,x,y=self.get_geometry(parent)
        widget.geometry('+%d+%d' %(x+200,y+250))
        return

class create_menubutton:

    def __init__(self,window,row=0,column=0,l=[],start_label=''):
        """Create a MenuButton from a list of [label,value,command] entries. """
        self.var=StringVar()
        self.var.set(start_label)
        self.button=Menubutton(window,textvariable=self.var,relief=RAISED)
        self.menu=Menu(self.button,tearoff=0)
        self.button['menu']=self.menu
        #
        # Add the labels and the commands
        #
        for label,value,command in l:
            self.menu.add_radiobutton(label=label,
                                      variable=self.var,
                                      value=value,
                                      indicatoron=0,
                                      command=command)
        self.button.grid(row=row,column=column,sticky='news')
        return

    #
    # -----
    #


class progress_window:

    def __init__(self,window,title):
        """Open a progress window"""
        self.parent_win=window
        #
        self.progress_win=Toplevel()
        self.progress_win.transient(self.parent_win)
        self.progress_win.focus_set()
        self.progress_win.title(title)
        self.progress_win.geometry('+300+150')
        self.prog_x=300
        self.prog_can=Canvas(self.progress_win,bg='white',width=self.prog_x,height=60)
        self.prog_can.grid(row=0,column=0)
        #
        self.progress=IntVar()
        self.xl=Label(self.progress_win,textvariable=self.progress,bg='white')
        self.xl.grid(row=1,column=0,sticky='news')
        self.box=None
        return

    def update_progress(self,fraction,level=1):

        y_pos=(level-1)*40+20
        bar_length=(self.prog_x-20)*fraction
        old_obj=self.box
        self.box=self.prog_can.create_polygon(10,y_pos,
                                              10+bar_length,y_pos,
                                              10+bar_length,y_pos+20,
                                              10,y_pos+20,fill='blue',
                                              outline='black')
        #
        # Text
        #
        txt='%4.1fx complete' %(100*fraction)
        import string
        txt=string.replace(txt,'x','%')
        self.progress.set(txt)
        if old_obj:
            self.prog_can.delete(old_obj)
        else:
            if level==1:
                bar_length=self.prog_x-20
                self.prog_can.create_polygon(10,y_pos,
                                             10+bar_length,y_pos,
                                             10+bar_length,y_pos+20,
                                             10,y_pos+20,fill='white',
                                             outline='black')
        self.progress_win.update()
        #self.parent_win.update()
        return

    def close(self):
        """Close the window"""
        self.progress_win.destroy()
        return

    def __del__(self):
        """Destructor"""
        self.close()
        return

