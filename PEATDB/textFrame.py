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
import Pmw

class textFrame:
    """Class to provide a frame for text editing"""
    def __init__(self,parent,filename=None,
                    title=None,
                    parent_info=None,
                    w=80,h=24):
        self.frame=Toplevel()
        self.frame.geometry('+400+400')
        if title==None:
            self.frame.title('Text Area')
        else:
            self.frame.title(title)
        if parent:
            self.parent=parent
        if parent_info:
            self.parent_info = parent_info
        row=0
        column=0
        details_width=w
        details_height=h
        import tkFont
        thefont = tkFont.Font (family="Arial", size=12)

        yscrollbar=Scrollbar(self.frame,orient='vertical',width=14)
        yscrollbar.grid(row=row,column=3,sticky='news',padx=2)
        xscrollbar=Scrollbar(self.frame,orient='horizontal',width=14)
        xscrollbar.grid(row=row+1,column=0,columnspan=3,sticky='news')
        self.details=Text(self.frame,background='white',
                      foreground='black',
                      state=NORMAL,
                      exportselection=1,
                      yscrollcommand=yscrollbar.set,
                      xscrollcommand=xscrollbar.set,
                      wrap='word',
                      font=thefont,
                      bg='lightgray')
        self.details.config(width=w)
        self.details.config(height=h)
        yscrollbar.config(command=self.details.yview)
        xscrollbar.config(command=self.details.xview)
        self.details.grid(row=0,column=0,columnspan=3,sticky='NEWS',padx=2,pady=2)
        self.details.config(state=NORMAL)
        self.frame.rowconfigure(0,weight=1)
        self.frame.columnconfigure(0,weight=1)
        #
        # Bottom panel consisting of close and save buttons
        #
        buttons=[['Save',self.save],
                 ['Save As',self.saveas],
                 ['Close',self.close]]
        row=2
        column=0
        for button,command in buttons:
            x=Button(self.frame,text=button,command=command)
            x.grid(row=row,column=column,sticky='NEWS',padx=3,pady=3)
            column=column+1

        self.frame.bind("<Destroy>", self.close )
        self.details.bind('<KeyPress>' , self.getcontents)
        return

    def load_from_file(self,filename):
        """Displays text from a simple text file"""
        if not filename:
            return
        import fileinput
        for line in fileinput.input(filename):
            self.details.insert(END, line)
        self.current_contents = self.details.get(1.0, END)
        self.filename = filename
        return


    def load_text(self,textdump):
        """Displays a text dump in the text widget, line by line"""
        self.details.delete(1.0, END)
        if textdump:
             for line in textdump:
                str = ''
                for S in line:
                    str = str+S
                self.details.insert(END, str)
        self.current_contents = self.details.get(1.0, END)
        return

    def load_text1(self,textdump,clear=0):
        """Displays a text dump in the text widget, line by line"""
        if clear==1:
            self.details.delete(1.0, END)
        if textdump:
            for line in textdump:
                str = ''
                for S in line:
                    str = str+S
                self.details.insert(END, (str+'\n'))
        self.current_contents = self.details.get(1.0, END)
        return

    def add_text(self,str):
        """Append some text item to the end of the list"""
        self.details.insert(END, '\n')
        self.details.insert(END, str)
        self.current_contents = self.details.get(1.0, END)
        return

    def add_text1(self,str):
        """Append some text item to the end of the list"""
        self.details.insert(END, '\n')
        self.details.insert(END, str)
        self.current_contents = self.details.get(1.0, END)
        return

    def load_list(self, inputlist):
        """Displays a list the text widget, element-wise"""
        self.details.delete(1.0, END)
        if inputlist:
             for line in inputlist:
                 for item in line:
                     s=str(item)
                     self.details.insert(END, s)
                     self.details.insert(END, ' ')
                 self.details.insert(END, '\n')
        self.current_contents = self.details.get(1.0, END)
        return

    def resize(self,event):
        """Make sure the visible portion of the text widget is resized"""
        if event.widget==self.frame:
            Y=event.height
            X=event.width
            self.details.configure(width=X-2, height=Y-2)
        return

    def doSave(self, filename):
        #self.getcontents()
        self.current_contents = self.details.get(1.0, END)
        try:
            fd=open(filename,'w')
            fd.write(self.current_contents)
            fd.close()
            print "File written ",filename
            self.filename = filename
        except:
            print "error: could not write file ",filename
        return    
            
    def save(self):
        import os
        if hasattr(self, 'filename') and self.filename!=None:
            if not os.path.exists(self.filename):
                return
            self.doSave(self.filename)
        else:
            self.saveas()
        return
    
    def saveas(self,filename=None):
        """Save the contents of the text widget to a file"""
        import tkFileDialog, os
        savedir=os.getcwd()
        if filename==None:
            filename=tkFileDialog.asksaveasfilename(parent=self.frame,
                                     defaultextension='.txt',
                                     initialdir=savedir,
                                     filetypes=[("Text files","*.txt"),("All files","*.*")])

        if not filename:
            return
        self.doSave(filename)
        return

    def getcontents(self,event=None):
        """Get the contents of the text frame"""
        self.current_contents = self.details.get(1.0, END)

        return

    def close(self,event=None):
        """Close the frame"""
        if hasattr(self,'parent_info'):
            self.return_peat_data()
        self.frame.destroy()
        #destroy parent instance
        #self.parent.text_frame=None
        return

    #method for returning data to peat
    def return_peat_data(self):
        try:
            self.getcontents()
        except:
            pass
        notesdata={}
        notesdata['parent_info']=self.parent_info
        notesdata['text']=self.current_contents
        self.parent.notesdata=notesdata.copy()
        return

class HyperlinkManager:

    def __init__(self, text):

        self.text = text

        self.text.tag_config("hyper", foreground="blue", underline=1)

        self.text.tag_bind("hyper", "<Enter>", self._enter)
        self.text.tag_bind("hyper", "<Leave>", self._leave)
        self.text.tag_bind("hyper", "<Button-1>", self._click)

        self.reset()

    def reset(self):
        self.links = {}

    def add(self, action):
        # add an action to the manager.  returns tags to use in
        # associated text widget
        tag = "hyper-%d" % len(self.links)
        self.links[tag] = action
        return "hyper", tag

    def add_link(self, link):
        # add an action to the manager.  returns tags to use in
        # associated text widget
        tag = "hyper-%d" % len(self.links)
        def action():
            import webbrowser
            webbrowser.open(link,autoraise=1)
        self.links[tag] = action
        return "hyper", tag

    def _enter(self, event):
        self.text.config(cursor="hand2")

    def _leave(self, event):
        self.text.config(cursor="")

    def _click(self, event):
        for tag in self.text.tag_names(CURRENT):
            if tag[:6] == "hyper-":
                self.links[tag]()
                return

