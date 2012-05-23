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

class DNA_Edit:

    def find_seq(self,event=None):
        #
        # Open a find box
        #
        self.find_win=Toplevel()
        self.find_win.title('Find DNA sequence')
        self.find_win.geometry('+%d+%d' %(self.master.winfo_rootx()+50,
                                          self.master.winfo_rooty()+50))
        #
        self.find_var=StringVar()
        lab=Label(self.find_win,text='Find:')
        lab.grid(row=0,column=0,sticky='news')
        Ebox=Entry(self.find_win,textvariable=self.find_var,width=20)
        Ebox.grid(row=0,column=1,sticky='news')
        #
        #
        findnext=Button(self.find_win,text='Find next',command=self.find_next)
        findnext.grid(row=1,column=0,sticky='news')
        close=Button(self.find_win,text='Close',command=self.close_find)
        close.grid(row=1,column=1,sticky='news')
        self.find_status=StringVar()
        self.find_status.set('No string entered')
        stat=Label(self.find_win,textvariable=self.find_status)
        stat.grid(row=2,column=0,columnspan=2,sticky='news')
        #
        #
        self.find_index=None
        #
        #
        self.find_win.bind('<KeyPress>',self.update_find)
        self.find_win.grab_set()
        Ebox.focus_set()
        return

    #
    # ---------
    #

    def update_find(self,event):
        #
        # find the string
        #
        ok_target=self.get_ok_DNAseq()
        if ok_target=='':
            return
        import string
        pos=string.find(self.data['DNAseq'],string.upper(ok_target))
        if pos==-1:
            self.find_index=None
            self.find_status.set('Not found')
            #
            # Clear any marked sequence
            #
            if self.data.has_key('sel_objs'):
                for obj in self.data['sel_objs'].keys():
                    self.seqframe.delete(obj)
            return
        self.find_status.set('Found, position: %s' %pos)
        self.mark_base_selection(pos,pos+len(ok_target))
        self.find_index=pos
        self.center_on_find(pos)
        return


    def get_ok_DNAseq(self):
        #
        # Clean the sequence entered
        #
        target=self.find_var.get()
        DNA=['A','C','G','T','a','c','g','t']
        ok_target=''
        for letter in target:
            if letter in DNA:
                ok_target=ok_target+letter
        self.find_var.set(ok_target)
        return ok_target


    def find_next(self,event=None):
        #
        # Find the next occurence
        #
        if self.find_index!=None:
            old_find=self.find_index
            ok_target=self.get_ok_DNAseq()
            import string
            nextpos=self.data['DNAseq'][old_find+1:].find(ok_target.upper())
            if nextpos==-1:
                nextpos=string.find(self.data['DNAseq'],string.upper(ok_target))-1
                old_find=0
                #self.find_status.set('No more occurences')
                #return
            #
            # Calc the real pos
            #
            pos=nextpos+old_find+1
            self.find_status.set('Found, position: %s' %pos)
            self.mark_base_selection(pos,pos+len(ok_target))
            self.find_index=pos
            self.center_on_find(pos)
        return

    #
    # --------
    #

    def close_find(self,event=None):
        #
        # Clear any marked sequence
        #
        if self.data.has_key('sel_objs'):
            for obj in self.data['sel_objs'].keys():
                self.seqframe.delete(obj)
        #
        # Close the find window
        #
        self.find_win.unbind('<KeyPress>')
        self.find_win.destroy()
        return

    #
    # -------
    #

    def center_on_find(self,position):
        x,y=self.get_base_pos_on_screen(position)
        #
        # We need to center on the base
        #
        center=max(0.0,x-(self.x_size-self.canvas_border_x)/2.0)
        #
        # Moveto works in fractions of the screen, so get the fraction
        #
        frac=center/self.canvas_x
        self.seqframe.xview('moveto', frac)
        self.seqframe.yview('moveto', 0.2)
        return


    def start_selection(self,event=None):
        """Start selecting a sequence"""

        if not self.data['DNAseq']:
            return
        self.clear_selection()
        # find the base clicked
        x=self.seqframe.canvasx(event.x)
        y=self.seqframe.canvasy(event.y)
        position=self.get_DNApos_fromcoords(x,y)
        if position==None or position<0 or position>len(self.data['DNAseq']):
            return

        # Start the selection
        self.data['DNA_selection']={}
        self.data['sel_objs']={}

        self.data['DNA_selection']['start']=position
        self.data['DNA_selection']['stop']=position
        self.data['DNA_selection']['active']=1

        # Mark it
        self.mark_base_selection(self.data['DNA_selection']['start'],self.data['DNA_selection']['stop'])
        self.seqframe.bind('<B1-Motion>',self.extend_selection)
        #self.seqframe.bind('<Shift-Button-1>',self.extend_selection)

        # Initialise counter
        self.DNAfragment_counter.set(abs(self.data['DNA_selection']['stop']-self.data['DNA_selection']['start']))

        if getattr(self,'routine_call_DNAselection',None):
            (self.routine_call_DNAselection)()
        return


    def end_selection(self,event=None):
        """End selecting sequence"""

        if not self.data['DNAseq']:
            return

        # End selection
        if self.data.has_key('DNA_selection'):
            if self.data['DNA_selection']['active']:
                self.data['DNA_selection']['active']=None
                self.seqframe.unbind("<B1-Motion>")
                #self.seqframe.unbind('<Shift-Button-1>')

        # Activate functions
        states={'Cut':1,
                'Copy':1,
                'Paste':1,
                'Delete':1}

        # Don't activate Paste unless we have something on the clipboard
        try:
            x=self.master.selection_get(selection = "CLIPBOARD")
        except TclError:
            states['Paste']=None
        self.set_pulldown_and_shortcut_states(states,self.seqmenu)
        return


    def extend_selection(self,event=None):
        """Extend the selection"""

        if self.data.has_key('DNA_selection'):
            if self.data['DNA_selection']['active']:
                x=self.seqframe.canvasx(event.x)
                y=self.seqframe.canvasy(event.y)
                position=self.get_DNApos_fromcoords(x,y)
                if position>=0 and position<len(self.data['DNAseq']):
                    self.data['DNA_selection']['stop']=position
                    self.mark_base_selection(self.data['DNA_selection']['start'],self.data['DNA_selection']['stop'])

                # Set the counter
                self.DNAfragment_counter.set(abs(self.data['DNA_selection']['stop']-self.data['DNA_selection']['start']))

                # Should we scroll?
                if event.x>self.x_size-self.canvas_border_x/2:
                    self.seqframe.xview('scroll',1,'units')
                elif event.x<5:
                    self.seqframe.xview('scroll',-11,'units')

        # Does anyone else want to know about this?
        if getattr(self,'routine_call_DNAselection',None):
            (self.routine_call_DNAselection)()
        return

    def clear_selection(self):
        """Clear DNA seq selection"""
        if self.data.has_key('DNA_selection'):
            del self.data['DNA_selection']

        # Clear display
        if self.data.has_key('sel_objs'):
            for obj in self.data['sel_objs'].keys():
                self.seqframe.delete(obj)
        return

    def mark_base_selection(self,start,stop,freeze=None):
        """Mark base selection"""

        if not self.data.has_key('sel_objs'):
            self.data['sel_objs']={}

        for obj in self.data['sel_objs'].keys():
            self.seqframe.delete(obj)

        # Draw the box
        bstart=min(start,stop)
        bend=max(start,stop)
        startx,y=self.get_base_pos_on_screen(bstart+1)
        endx,y=self.get_base_pos_on_screen(bend+1)
        colour='yellow'
        if start==stop:
            endx=endx+1
            startx=startx-1
            colour='black'
        ysize=5
        obj=self.seqframe.create_polygon(startx,y-ysize,
                                         endx,y-ysize,
                                         endx,y+ysize,
                                         startx,y+ysize,fill=colour)
        self.data['sel_objs'][obj]=1
        if freeze:
            obj=self.seqframe.create_polygon(startx,y-ysize,
                                             endx,y-ysize,
                                             endx,y+ysize,
                                             startx,y+ysize,
                                             outline='black',fill="")
            self.data['sel_objs'][obj]=1

        # Mark the bases that are selected
        font = self.getCurrentFont()
        for position in range(bstart,bend):
            markx,marky=self.get_base_pos_on_screen(position+1)
            obj=self.seqframe.create_text(markx,marky,text=self.data['DNAseq'][position],fill='red',font=font,anchor='w')
            self.data['sel_objs'][obj]=1
        return


    def cut_DNA(self,event=None,delete=None):
        """Cut the selected DNA stretch"""

        if self.data.has_key('DNA_selection'):
            if not self.data['DNA_selection']['active']:
                start=self.data['DNA_selection']['start']
                stop=self.data['DNA_selection']['stop']
                c_start=min(start,stop)
                c_end=max(start,stop)
                if not delete:
                    self.data['DNAclipboard']=self.data['DNAseq'][c_start:c_end]
                self.data['DNAseq']=self.data['DNAseq'][:c_start]+self.data['DNAseq'][c_end:]

                # Add to clipboard
                if not delete:
                    self.master.clipboard_clear()
                    self.master.clipboard_append(self.data['DNAclipboard'])

                # Clear the selection
                self.clear_selection()
                # Update everything
                self.update_sequence_window()
        return


    def delete_DNA(self,event=None):
        """Call a slightly different Cut function"""
        self.cut_DNA(delete=1)
        return

    def copy_DNA(self,event=None):
        #
        # Cut the selected DNA stretch
        #
        if self.data.has_key('DNA_selection'):
            if not self.data['DNA_selection']['active']:
                start=self.data['DNA_selection']['start']
                stop=self.data['DNA_selection']['stop']
                c_start=min(start,stop)
                c_end=max(start,stop)
                self.data['DNAclipboard']=self.data['DNAseq'][c_start:c_end]
                #
                # Mark it
                #
                self.mark_base_selection(self.data['DNA_selection']['start'],
                                         self.data['DNA_selection']['stop'],freeze=1)
                #
                # Add to clipboard
                #
                self.master.clipboard_clear()
                self.master.clipboard_append(self.data['DNAclipboard'])
        return


    def paste_DNA(self,event=None):
        """Paste stuff"""

        text=self.master.selection_get(selection = "CLIPBOARD")

        # Check that we have DNA

        import string
        text=string.replace(text,'\n','')
        text=string.replace(text,'\r','')
        text=string.replace(text,' ','')
        DNA=['A','C','G','T','a','c','g','t']
        for letter in text:
            if not letter in DNA:
                import tkMessageBox
                tkMessageBox.showwarning('Invalid DNA sequence',
                                         '"%s" is not a valid DNA sequence. \n"%s" is first illegal base.' %(text,letter),
                                         parent=self.master)
                return

        # Make upper case
        text=string.upper(text)

        # ok

        start=self.data['DNA_selection']['start']
        stop=self.data['DNA_selection']['stop']
        p_start=min(start,stop)
        p_end=max(start,stop)
        self.data['DNAseq']=self.data['DNAseq'][:p_start]+text+self.data['DNAseq'][p_end:]
        self.clear_selection()
        self.update_sequence_window()
        return


