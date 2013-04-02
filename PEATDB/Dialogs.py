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
import tkSimpleDialog, tkFileDialog, tkMessageBox

class myDialogs:

    def askopenfilename(self,defaultextension=None,initialdir=None,filetypes=None):
        #
        # Ask for the open filename
        #
        dir_frompref=None
        if not initialdir:
            initialdir=self.preferences.get('datadir')
            dir_frompref=1
        #
        # Call the Tkinter routine
        #
        import tkFileDialog
        filename=tkFileDialog.askopenfilename(defaultextension=defaultextension,
                                              initialdir=initialdir,
                                              filetypes=filetypes)
        #
        # Check if the user changed dir
        #
        import os
        if filename and dir_frompref:
            ldir=os.path.split(filename)[0]
            if ldir!=initialdir:
                self.preferences.set('datadir',dir)
        #
        # Return
        #
        return filename


class entry_field:
    """ create a simple label - entry field pair"""

    def __init__(self,window,row,column,name,ftype='text',items=[],default_item=''):
        """Create a simple entry field in the window"""
        self.label=Label(window,text=name)
        self.label.grid(row=row,column=column,sticky='news')
        #
        # Now the entry itself
        #
        if ftype=='text' or ftype=='number':
            self.var=StringVar()
            self.entry=Entry(window,width=15,bg='white',textvariable=self.var)
            self.entry.grid(row=row,column=column+1,sticky='news')
        elif ftype=='menu':
            #
            # Pull-down like menu button
            #
            self.var=StringVar()
            self.var.set(default_item)
            self.button=Menubutton(window,textvariable=self.var,relief=RAISED)
            self.menu=Menu(self.button,tearoff=0)
            self.button['menu']=self.menu
            #
            # Add the options
            #
            for option in items:
                self.menu.add_radiobutton(label=option,
                                          variable=self.var,
                                          value=option,
                                          indicatoron=1)
            self.button.grid(row=row,column=column+1,sticky='news')
        return

class askyesnocancel(Frame):

    def __init__(self,title='Dialog',message='Message',parent=None):
        """Display a Yes,No,Cancel dialog"""
        self.parent=parent
        self.win=Toplevel()
        self.win.transient(parent)
        self.win.title(title)
        self.add_buttons(message)
        self.win.grab_set()
        self.win.protocol('WM_DELETE_WINDOW',self.cancel)
        self.win.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                      parent.winfo_rooty()+50))


        #self.win.initial_focus.focus_set()
        self.win.bind("&lt;Return>", self.yes)
        self.win.bind("&lt;Escape>", self.cancel)
        self.win.wait_window(self.win)
        return None

    def add_buttons(self,message):
        """Add the buttons"""
        Label(self.win,text=message).pack(side=TOP,fill=X)
        Button(self.win,text='Yes',command=self.yes).pack(side=LEFT,fill=X,expand=1)
        Button(self.win,text='No',command=self.no).pack(side=LEFT,fill=X,expand=1)
        Button(self.win,text='Cancel',command=self.cancel).pack(side=LEFT,fill=X,expand=1)
        return

    def yes(self):
        """Return teh yes value"""
        self.closewin()
        self.result='yes'
        return

    def no(self):
        """Return the no value"""
        self.closewin()
        self.result='no'
        return

    def cancel(self):
        """Return none"""
        self.result=None
        self.parent.focus_set()
        self.win.destroy()
        return

    def closewin(self,event=None):
        """Close the window"""
        self.win.destroy()
        self.win.update_idletasks()
        return

class PEATDialog(Tk):

    def __init__(self, parent, option='generic',message='Working'):
        if parent!=None:
            self.main=parent
            self.rootx=self.main.winfo_rootx()
            self.rooty=self.main.winfo_rooty()
            if option=='generic':
                self.show_loading_progress(message)
            elif option=='pmw':
                self.show_pmw_dialog()
            elif option == 'progressbar':
                self.show_progressbar(message)
            elif option == 'busybar':
                self.show_busy_bar(message)
            self.themessage=message
        else:
            return

    def show_loading_progress(self, message='Working..'):
        """Show plain progress window for loading of data"""
        self.progress_win=Toplevel() # Open a new window
        self.progress_win.title("Please Wait")
        self.progress_win.geometry('+%d+%d' %(self.rootx+200,self.rooty+200))
        #force on top
        self.progress_win.grab_set()
        self.progress_win.transient(self.main)
        progrlbl = Label(self.progress_win,text=message)
        progrlbl.grid(row=0,column=0,sticky='news',padx=6,pady=4)
        progrlbl2 = Label(self.progress_win,text='Please Wait..')
        progrlbl2.grid(row=1,column=0,sticky='news',padx=6,pady=4)
        self.progress_win.update()
        return

    def show_progressbar(self,message):
        """Show progress bar window for loading of data"""
        self.progress_win=Toplevel() # Open a new window
        self.progress_win.title("Please Wait")
        self.progress_win.geometry('+%d+%d' %(self.rootx+200,self.rooty+200))
        #force on top
        try:
            self.progress_win.grab_set()
            self.progress_win.transient(self.main)
        except:
            pass
        lbl = Label(self.progress_win,text=message,font='Arial 16')
        #lbl = Label(self.progress_win,text=self.themessage,font='Arial 16')
        lbl.grid(row=0,column=0,columnspan=2,sticky='news',padx=6,pady=4)
        progrlbl = Label(self.progress_win,text='Progress:')
        progrlbl.grid(row=1,column=0,sticky='news',padx=2,pady=4)
        import ProgressBar
        self.bar = ProgressBar.ProgressBar(self.progress_win)
        self.bar.frame.grid(row=1,column=1,columnspan=2,padx=2,pady=4)
        return

    def show_busy_bar(self,message='Please Wait'):
        """Show busy bar window for indeterminate time"""
        self.progress_win=Toplevel() # Open a new window
        self.progress_win.title(message)
        self.progress_win.geometry('+%d+%d' %(self.rootx+200,self.rooty+200))
        #force on top
        self.progress_win.grab_set()
        self.progress_win.transient(self.main)
        import BusyBar
        l = Label(self.progress_win, text=message)
        l.grid(row=0,column=0,sticky='nw',padx=6,pady=6)
        self.bb = BusyBar.BusyBar(self.progress_win, text='Please Wait')
        self.bb.grid(row=1,column=0,sticky='nw',padx=6,pady=6)
        return

    def update_busy(self):
        self.bb.update()
        return

    def busy_on(self):
        if hasattr(self,'bb'):
            self.bb.on()
        return

    def busy_off(self):
        if hasattr(self,'bb'):
            self.bb.of()
        return

    def update_progress(self,value):
        """update the progress bar if correct type"""
        if hasattr(self,'bar'):
            self.bar.updateProgress(newValue=value)
            self.progress_win.update_idletasks()
            #self.master.update_idletasks()
        else:
            return
        return

    @classmethod
    def createBusyBar(cls, parent):
        import BusyBar
        bb = BusyBar.BusyBar(parent, text='Busy')
        return bb

    @classmethod
    def createProgressBar(cls, parent):
        import ProgressBar
        bar = ProgressBar.ProgressBar(parent)
        return bar

    def show_pmw_dialog(self):
        """Show a pmw dialog window"""
        import Pmw
        # Create the dialog.
        self.dialog = Pmw.Dialog(self.main,
            buttons = ('OK', 'Apply', 'Cancel', 'Help'),
            defaultbutton = 'OK',
            title = 'My dialog')
            #command = self.execute)
        self.dialog.withdraw()

        # Add some contents to the dialog.
        w = Label(self.dialog.interior(),
            #text = self.themessage,
            text='Working',
            background = 'black',
            foreground = 'white',
            pady = 20)
        w.pack(expand = 1, fill = 'both', padx = 4, pady = 4)

        # Create the window excluded from showbusycursor.
        self.excluded = Pmw.MessageDialog(self.main,
            title = 'I still work',
            message_text =
                'This window will not get\n' +
                'a busy cursor when modal dialogs\n' +
                'are activated.  In addition,\n' +
                'you can still interact with\n' +
                'this window when a "no grab"\n' +
                'modal dialog is displayed.')
        self.excluded.withdraw()
        Pmw.setbusycursorattributes(self.excluded.component('hull'),
            exclude = 1)
        return


    def close(self):
        if self.progress_win:
            self.progress_win.destroy()
        elif self.dialog:
            self.dialog.destroy()
        return

def createBusywin(parent, message):
    """Show busy bar window for indeterminate time"""
    rootx = parent.winfo_rootx()
    rooty = parent.winfo_rooty()
    progress_win=Toplevel()
    progress_win.title("Please Wait")
    progress_win.geometry('+%d+%d' %(rootx+200,rooty+200))
    progress_win.grab_set()
    progress_win.transient(parent)
    import BusyBar
    l = Label(progress_win, text=message)
    l.grid(row=0,column=0,sticky='nw',padx=6,pady=6)
    bb = BusyBar.BusyBar(progress_win, text='Please Wait')
    bb.grid(row=1,column=0,sticky='nw',padx=6,pady=6)
    return progress_win, bb

class TopLevelModalDialog(Toplevel):
    def __init__(self, parent, width=300, height=300):
        Toplevel.__init__(self, parent)
        self.transient(parent)
        self.title('PEATSA Results Export')
        self.parent = parent
        self.geometry('%sx%s+300+200' %(width,height))
        self.body = Frame(self)
        self.initial_focus = self.body
        self.body.pack(fill=BOTH,expand=1,padx=5, pady=5)
        self.grab_set()
        if not self.initial_focus:
            self.initial_focus = self
        self.protocol("WM_DELETE_WINDOW", self.cancel)
        return

    def cancel(self):
        self.parent.focus_set()
        self.destroy()
        return

class MultipleValDialog(tkSimpleDialog.Dialog):
    """Simple dialog to get multiple values"""

    def __init__(self, parent, title=None, initialvalues=None, labels=None, types=None):
        if labels != None and types != NoneType:
            self.initialvalues = initialvalues
            self.labels = labels
            self.types = types
        tkSimpleDialog.Dialog.__init__(self, parent, title)

    def body(self, master):

        r=0
        self.vrs=[];self.entries=[]
        for i in range(len(self.labels)):
            Label(master, text=self.labels[i]).grid(row=r, column=0,sticky='news')
            if self.types[i] == 'int':
                self.vrs.append(IntVar())
            else:
                self.vrs.append(StringVar())
            if self.types[i] == 'password':
                s='*'
            else:
                s=None

            if self.types[i] == 'list':
                button=Menubutton(master, textvariable=self.vrs[i],relief=RAISED)
                menu=Menu(button,tearoff=0)
                button['menu']=menu
                choices=self.initialvalues[i]
                for c in choices:
                    menu.add_radiobutton(label=c,
                                        variable=self.vrs[i],
                                        value=c,
                                        indicatoron=1)
                self.entries.append(button)
                self.vrs[i].set(self.initialvalues[i][0])
            else:
                self.vrs[i].set(self.initialvalues[i])
                self.entries.append(Entry(master, textvariable=self.vrs[i], show=s, bg='white'))
            self.entries[i].grid(row=r, column=1,padx=2,pady=2,sticky='news')
            r+=1

        return self.entries[0] # initial focus

    def apply(self):
        self.result = True
        self.results = []
        for i in range(len(self.labels)):
            self.results.append(self.vrs[i].get())
        return
