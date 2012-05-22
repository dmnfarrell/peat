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


"""DNAtool application
Authors: Jens Nielsen & Damien Farrell"""

import sys, os
import tkFont
import platform

# Append the path to PEAT_DB
DNAtooldir = os.path.split(__file__)[0]

#split1=os.path.split(sys.path[0])[0]
#sys.path.append(split1)
#sys.path.append(os.path.split(split1)[0])
#the above leads to opaque imports later in the code...bi

from Tkinter import *
import Pmw
from DNAtool_IO import *
from DNA_Edit import *
from evaluate_primer import *
from Restriction_Digest_driver import *
from primer_design_GUI import *
from PEATDB.GUI_helper import GUI_help
from PEATDB.tooltip import *
from PEATDB.Prefs import Preferences

class MainWindow(Frame, DNA_IO, DNA_Edit, Restriction_Digest,
                    primer_design_GUI, GUI_help, evaluate_primer):
    """Class providing main window - should be considered legacy code"""

    def __init__(self, parent=None, data_passed=None, openwin=1):

        """If we are main then get variables from PEAT"""
        self.currplatform=platform.system()
        self.parent=parent

        # Define the look and feel
        self.font="Times 12 bold"
        self.bg_colour="grey"
        self.fg_colour="black"
        self.x_size=1000
        self.y_size=400
        self.canvas_border_x=50
        self.canvas_border_y=140
        self.canvas_x=self.x_size-self.canvas_border_x
        self.canvas_y=self.y_size-self.canvas_border_y
        self.canvas_height=700
        self.seq_row=self.canvas_height/2

        # Positioning bases and aas on the screen
        self.seq_xstart=70
        #scale for spacing between bases
        self.base_scale=10
        #offset value of ORF for vertical spacing
        self.orf_offset=8
        #vertical scale value, to adjust for font sizes
        self.y_scale=1

        self.seqfont='Courier 10'
        self.labelfont='Courier 10'
        self.restr_font='Arial 10'
        self.backgrcolor='#f6f9f6'

        # Graphics objects
        self.seq_win_objs={}
        self.digest_win_objs={}
        #restriction site objects
        self.overview_win=None
        self.overview_win_objs={}

        # Windows
        self.digest_window=None

        # Data that we need
        self.data={}
        self.data['dnaseqfile']=None
        self.data['DNAseq']=None
        self.data['DNAseq_status']=None
        self.data['cut_pos']=None
        self.data['used_enzymes']=None
        self.data['ORFS']=None

        # Dir for Data
        import os
        self.data['datadir']=os.getcwd()

        # Other variables
        self.data['Project saved']=None
        self.data['Project filename']=None

        # Vars for the primer design
        self.detailed_objs={}
        self.pDB_open=None
        # Var for sequencing
        self.sequ_win=None
        self.maxseqlevel=0  #level of highest sequence
        self.primer_displayed=0

        # Get the data that we got from our parent
        if data_passed:
            import copy
            for key in data_passed.keys():
                self.data[key]=copy.deepcopy(data_passed[key])

        self.check_primers()        
        
        if openwin:
            self.do_window()            
            #load display preferences if they are there
            self.load_preferences()             
            self.main_pulldown()
            self.init_tkvars()
            self.init_bindings()

        # Open the window
        if openwin:
            # Update the view
            self.update_sequence_window()
            self.open_pDB()
            self.assess_status()
       
        return

    def check_primers(self):
        """Check that all primers are ok"""
        primer_keys=['sequence','description']
        if self.data.has_key('primer_dict'):
            for primer in self.data['primer_dict'].keys():
                for key in primer_keys:
                    if not self.data['primer_dict'][primer].has_key(key):
                        self.data['primer_dict'][primer][key]='None'

    def do_window(self):
        """Create the main window """
        if not self.parent:
            Frame.__init__(self)
            self.master_win = self.master
        else:
            self.master = Toplevel()
            self.master_win = self.parent.master
        print self.master, self.master_win
        
        self.master.title("DNA sequence manipulation")
        self.master.geometry("%dx%d+%d+%d" %(self.x_size,self.y_size,100,100))

        # Set up the main frame
        #self.master.rowconfigure(0,weight=1)
        #self.master.columnconfigure(0,weight=1)
        #self.master.grid(sticky=W)

        self.init_vars()

        #Text box for top
        label1=Label(self.master, text="File holding DNA sequence ",font=self.font)
        label1.grid(row=0,column=0, sticky=W)

        #Entry field
        self.filename=Entry(self.master)

        # If we have a filename then insert it
        if self.data['dnaseqfile']:
            self.filename.insert(INSERT,self.data['dnaseqfile'])
        self.filename.grid(row=0,column=1,sticky='we')


        # Button for selecting the pdb file
        self.loadbutton=Button(self.master,text='Browse',command=self.dnaseq_read,
                               font=self.font,fg=self.fg_colour,bg=self.bg_colour)
        self.loadbutton.grid(row=0,column=2,sticky='W')
        t1 = ToolTip(self.loadbutton, follow_mouse=1, text="Select pdb file")

        # Information on the PEAT_DB record connected with the data
        Label(self.master,text='PEAT_DB record:').grid(row=0,column=3,sticky=E)
        self.PEAT_DBrecord=StringVar()
        Label(self.master,textvariable=self.PEAT_DBrecord).grid(row=0,column=4)

        # Exit button
        exit_text='Exit'
        if self.parent:
            exit_text='Return to DataBase'
        exitbutton=Button(self.master,text=exit_text,command=self.quit)
        exitbutton.grid(row=1,column=0,sticky='WE')
        t2 = ToolTip(exitbutton, follow_mouse=1, text="Exits and closes window")

        # Button for selecting the ORF
        self.ORFbutton=Button(self.master,text='Select ORF',command=self.select_ORF)
        self.ORFbutton.grid(row=1,column=1,columnspan=1,sticky='WE')
        self.ORFbutton.configure(state=DISABLED)
        t3 = ToolTip(self.ORFbutton, follow_mouse=1, text="Select open reading frame")

        # Open primer Database
        self.pDB_but=Button(self.master,text="Primer Database",command=self.open_pDB)
        self.pDB_but.grid(row=1,column=2,sticky='news')

        # Button for opening overview window
        self.overview_button=Checkbutton(self.master,text='ORF Overview',command=self.overview_on_off)
        self.overview_button.grid(row=1,column=3,sticky='E')

        # Restriction digest overview
        self.restr_detail=IntVar()
        self.restr_button=Checkbutton(self.master,text='Digest details',
                                      command=self.restr_on_off,
                                      var=self.restr_detail,onvalue=1,offvalue=0)
        self.restr_button.grid(row=2,column=3,sticky='E')
        self.restr_detail.set(0)


        # Create the main sequence window - this is where we will have the DNA sequence,
        # the protein sequence, restriction sites etc

        lbl1=Label(self.master,text='Sequence window')
        lbl1.grid(row=2,column=0,sticky='W')

        # Counter for selecting DNA fragments
        Label(self.master,text='Size of selected DNA fragment',bg='yellow').grid(row=2,column=1,columnspan=1)
        self.DNAfragment_counter=IntVar()
        self.DNAfragment_counter.set(0)
        Label(self.master,textvariable=self.DNAfragment_counter,bg='yellow').grid(row=2,column=2,columnspan=1,sticky='W')

        # Scrollbars
        span=6
        scrollbar=Scrollbar(self.master,orient='horizontal')
        scrollbar.grid(row=4,column=0,columnspan=span,sticky='NEWS')
        yscrollbar=Scrollbar(self.master,orient='vertical')
        yscrollbar.grid(row=3,column=6,rowspan=3,sticky='NEWS')


        # Canvas, draws the sequence
        self.seqframe=Canvas(self.master,bd=8,bg='white',width=self.x_size-10, height=self.canvas_height,
                             xscrollcommand=scrollbar.set, yscrollcommand=yscrollbar.set,
                             scrollregion=(0,0,self.canvas_x,self.canvas_y),
                             xscrollincrement=0.0)
        self.seqframe.grid(row=3,column=0,columnspan=span)

        # Resizing
        self.master.bind("<Configure>",self.resize)

        # Destroying the window
        self.master.protocol("WM_DELETE_WINDOW",self.quit)

        # Activate scrollbar
        scrollbar.config(command=self.seqframe.xview)
        yscrollbar.config(command=self.seqframe.yview)
        self.seqframe.xview("moveto", 0)
        self.seqframe.yview("moveto", 0.2)
        self.seqframe.bind('<Button-4>', lambda event: event.widget.yview_scroll(-1, UNITS))
        self.seqframe.bind('<Button-5>', lambda event: event.widget.yview_scroll(1, UNITS))

        #self.master.geometry("%dx%d" %(self.x_size,self.y_size))
        #self.update_sequence_window()

        return

    def init_tkvars(self):
        """Initialise the Tk variables"""

        self.base_scale_input=IntVar()
        self.base_scale_input.set(10)
        self.seqfont_input=StringVar()
        self.seqfont_input.set('Courier')
        self.seqfontsize_input=IntVar()
        if 'Darwin' in self.currplatform:
            self.seqfontsize_input.set(16)
        else:
            self.seqfontsize_input.set(14)
        self.fontstyle_input=IntVar()
        self.fontstyle_input.set(1)
        self.restr_site_font_input=StringVar()
        self.restr_site_font_input.set('Arial')
        self.resnum=IntVar()
        self.resnum.set(1)
        # Method for calculating Tm of primers
        self.Tm_method=StringVar()
        self.Tm_method.set('Stratagene')
        return


    def resize(self,event):
        """
        Make sure the visible portion of the canvas is resized with the
        window
        """
        if event.widget==self.master:
            Y=event.height
            X=event.width
            self.seqframe.configure(width=X-self.canvas_border_x,
                              height=Y-self.canvas_border_y)
        return


    def main_pulldown(self):
        """
        Create the main pulldown menu
        """
        self.menu=Menu(self.master)
        # File menu
        #
        self.file_menu=Menu(self.menu,tearoff=0)

        #self.file_menu.add_separator()
        self.file_menu.add_command(label='Open Project',command=self.project_open)
        self.file_menu.add_command(label='Close Project',command=self.project_close)
        self.file_menu.add_command(label='Save Project',command=self.project_save)
        self.file_menu.add_command(label='Save Project As',command=self.project_saveas)
        self.file_menu.add_command(label='Exit',command=self.quit)
        self.menu.add_cascade(label='File',menu=self.file_menu)
        # -----------------------------------------------
        # Sequence menu
        #
        self.seqmenu={'01Open':{'cmd':self.dnaseq_read},
                      '02Save':{'cmd':self.save_dnaseq},
                      '03Save As':{'cmd':self.save_dnaseq},
                      '04sep':{None:None},
                      '05Find':{'cmd':self.find_seq,'sc':'Ctrl+F'},
                      '06Cut':{'cmd':self.cut_DNA,'sc':'Ctrl+X'},
                      '07Copy':{'cmd':self.copy_DNA,'sc':'Ctrl+C'},
                      '08Paste':{'cmd':self.paste_DNA,'sc':'Ctrl+V'},
                      '09Delete':{'cmd':self.delete_DNA,'sc':'Del'}}
        self.seqmenu=self.create_pulldown(self.menu,self.seqmenu)
        #
        # Invert sequence
        #
        self.invert_seq_var=IntVar()
        self.invert_seq_var.set(0)
        self.seqmenu['var'].add_checkbutton(label='Invert',
                                            command=self.invert_seq,
                                            variable=self.invert_seq_var,onvalue=1,offvalue=0)
        #
        # Complementary seq
        #
        self.complement_seq_var=IntVar()
        self.complement_seq_var.set(0)
        self.seqmenu['var'].add_checkbutton(label='Complementary',
                                            command=self.complementary_seq,
                                            variable=self.complement_seq_var,onvalue=1,offvalue=0)
        #
        # Set restriction digest parameters
        #
        self.seqmenu['var'].add_command(label="Configure restriction digest",command=self.win_select_enzymes)
        self.seqmenu['var'].add_separator()
        # Option to set spacing of bases, if needed
        self.seqmenu['var'].add_command(label="Sequence Display Setup",command=self.seq_display_settings)
        self.seqmenu['var'].add_separator()

        # Colour-coding of bases
        self.colour_seq_var=IntVar()
        self.colour_seq_var.set(2)

        # self.seqmenu['var'].add_checkbutton(label='Colour sequence',
                                            # command=self.update_sequence_window,
                                            # variable=self.colour_seq_var,onvalue=1,offvalue=0)

        self.seqmenu['var'].add_radiobutton(label="No Colour", command=self.update_sequence_window,
                                            variable=self.colour_seq_var, value=0)
        self.seqmenu['var'].add_radiobutton(label="Colour by base", command=self.update_sequence_window,
                                            variable=self.colour_seq_var, value=1)
        self.seqmenu['var'].add_radiobutton(label="Colour in threes", command=self.update_sequence_window,
                                            variable=self.colour_seq_var, value=2)

        #
        # Add the self.sequence_menu to the main menu
        #
        self.menu.add_cascade(label='DNA Sequence',menu=self.seqmenu['var'])
        # -----------------------------------------------
        # Primer menu
        #
        self.primer_menu=Menu(self.menu,tearoff=0)
        self.primer_menu.add_command(label='PCR primers',command=self.design_PCR_primer)
        self.primer_menu.add_command(label='Mutagenic primer',command=self.design_mutagenic_primer)
        self.primer_menu.add_command(label='Evaluate primer',command=self.do_evaluate_primer)
        self.primer_menu.add_command(label='About primer design',command=self.about_primer_design)

        self.menu.add_cascade(label='Primer Design',menu=self.primer_menu)


        # Add menu for checking sequencing runs
        self.show_comp_sequence=IntVar()

        #self.sequencing_menu=Menu(self.menu,tearoff=0)
        #self.sequencing_menu.add_checkbutton(label='Sequencing',
        #                                     command=self.sequencing_window,
        #                                     variable=self.sequencing_window_var,
        #                                     onvalue=1,offvalue=0)
        self.show_comp_sequence.set(0)
        #self.menu.add_cascade(label='Window',menu=self.sequencing_menu)

        # Add menu for importing and analysis of sequences
        self.analyse_sequences_menu=Menu(self.menu,tearoff=0)
        self.analyse_sequences_menu.add_command(label='Load Single File',
                                                command=self.compare_DNA_seq)
        self.analyse_sequences_menu.add_command(label='Load Mutiple Files',
                                                command=self.compare_multiple_seqs)
        self.analyse_sequences_menu.add_command(label='Clear Sequence(s)',
                                                command=self.clear_comp_seqs)
        self.menu.add_cascade(label='Sequence Analysis',menu=self.analyse_sequences_menu)


        #About menu
        self.help_menu=Menu(self.menu,tearoff=0)
        self.help_menu.add_command(label='About DNAtool',command=self.about_DNAtool)
        self.help_menu.add_command(label='Online documentation',command=self.online_documentation)
        self.menu.add_cascade(label='Help',menu=self.help_menu)

        self.master.config(menu=self.menu)

        # Set standard active and disabled states
        self.assess_status()
        return



    def assess_status(self):
        """Figures out which buttons to activate/deactivate"""

        # Update the PEATDB_record information
        if self.data.has_key('PEAT_DBrecord'):
            self.PEAT_DBrecord.set(self.data['PEAT_DBrecord'])

        # Do we have an open project?
        if self.data['Project filename']:
            self.file_menu.entryconfigure(0,state=DISABLED) # Open project
            self.file_menu.entryconfigure(1,state=NORMAL) # Close project
            self.file_menu.entryconfigure(3,state=NORMAL) # Save project as

            # Is the current project saved?
            if not self.data['Project saved']:
                self.file_menu.entryconfigure(2,state=NORMAL) # Save Project
            else:
                self.file_menu.entryconfigure(2,state=DISABLED) # Save Project
        else:

            # No project
            self.file_menu.entryconfigure(1,state=DISABLED) # Close project
            self.file_menu.entryconfigure(2,state=DISABLED) # Save Project
            self.file_menu.entryconfigure(3,state=DISABLED) # Save project as

        # Do we have a DNA sequence?
        if self.data['DNAseq']:
            self.ORFbutton.configure(state=NORMAL)     # ORF button
            # Sequence menu
            states={'Open':None,
                    'Save':1,
                    'Save As':1,
                    'Complementary':1,
                    'Inverse':1,
                    'Find':1}
            self.set_pulldown_and_shortcut_states(states,self.seqmenu)
            # Primer menu
            self.primer_menu.entryconfigure(0,state=ACTIVE) # PCR primer
            # .
            self.file_menu.entryconfigure(0,state=DISABLED) # Open project
            self.file_menu.entryconfigure(3,state=NORMAL) # Save project as
        else:
            # No DNA
            states={'Open':1,
                    'Save':None,
                    'Save As':None,
                    'Complementary':None,
                    'Inverse':None,
                    'Find':None,
                    'Cut':None,
                    'Copy':None,
                    'Paste':None,
                    'Delete':None}
            self.set_pulldown_and_shortcut_states(states,self.seqmenu)

            # Primer menu
            self.primer_menu.entryconfigure(0,state=DISABLED) # PCR primers
            self.primer_menu.entryconfigure(1,state=DISABLED) # Mutagenic primer

        # Has an ORF been selected?
        if self.data.has_key('ORF_selected'):
            self.primer_menu.entryconfigure(1,state=ACTIVE) # Mutagenic primer

            start=self.data['ORF_selected']['start']
            stop=self.data['ORF_selected']['stop']
        return

    def init_bindings(self):
        self.seqframe.bind("<Button-1>",self.handle_left_click)
        self.seqframe.bind("<ButtonRelease-1>",self.handle_left_release)
        self.seqframe.bind("<Button-3>",self.handle_right_click)
        self.seqframe.bind("<ButtonRelease-3>",self.handle_right_release)
 
        self.seqframe.bind("<Shift-Button-1>", self.handle_left_shift_click)
        self.seqframe.bind("<B3-Motion>", self.handle_right_motion)

        #self.seqframe.bind("<Button-3>", self.show_item)
        #self.seqframe.bind("<ButtonRelease-3>",self.remove_label)
        #self.seqframe.bind_all("<Control-KeyPress-x>",self.cut_DNA)  # Cut
        #self.seqframe.bind_all("<Control-KeyPress-c>",self.copy_DNA) # Copy
        #self.seqframe.bind_all("<Control-KeyPress-v>",self.paste_DNA) # Paste
        #self.seqframe.bind_all("<Delete>",self.delete_DNA) # Delete
        return

    def handle_left_click(self,event):
        """handle left mouse press on canvas"""
        c = self.seqframe
        if 'textlabel' in c.gettags(CURRENT):
            self.show_item(event)
        elif 'comparison_seq' in c.gettags(CURRENT):
            self.show_sequence_label(event)
        else:
            self.start_selection(event)

    def handle_left_release(self,event):
        """handle left mouse release on canvas"""
        c = self.seqframe
        if 'textlabel' in c.gettags(CURRENT):
            self.remove_recog_label(event)
        elif 'comparison_seq' in c.gettags(CURRENT):
            self.remove_seq_label(event)
        #else:
        #    self.end_selection(event)
        return

    def handle_left_shift_click(self, event):
        """handle shift left click for seq selection"""
        self.extend_selection(event)
        #self.end_selection(event)
        return

    def handle_right_click(self, event):
        """mainly for restr label movement"""
        c=self.seqframe
        if 'textlabel' in c.gettags(CURRENT):
            self.currobjs = c.find_withtag(CURRENT)            
        return
    
    def handle_right_release(self, event):
        """mainly for restr label movement"""
        if not hasattr(self,'currobjs') or len(self.currobjs)==0:
            return
        c=self.seqframe       
        y=self.seqframe.canvasy(event.y)
        for item in self.currobjs:            
            self.move_restriction_label(item, y)
        self.currobjs=[]    
        return
    
    def handle_right_motion(self, event):
        if len(self.currobjs)==0:
            return
        x=self.seqframe.canvasy(event.x)
        y=self.seqframe.canvasy(event.y)
        
        return
        
    def show_sequence_label(self, event):
        """how sequence name if it's a seq label, called from DNAtool main"""
        c=self.seqframe
        box = c.bbox(CURRENT)
        x1=box[0]
        y1=box[1]
        x2=box[2]
        y2=box[3]
        items=[]
        #make selection rectangle one pixel larger to include rect and text
        items=c.find_enclosed(x1-1,y1-1,x2+1,y2+1)

        import tkFont
        sfont = tkFont.Font (family='Arial', size=12,weight='bold')
        for obj in items:
            c.tag_raise(obj)
            #if item is text, get recog sequence and display
            for name in c.gettags(obj):
                #ignore other tags, just get name of seq
                if name!='current' and name!='comparison_seq':
                    obj=c.create_text(x2+3,y1-3,text=name,tags='seqlabel',
                                    font=sfont,width=120,anchor='nw')
                    box = c.bbox(obj)
                    rect = c.create_rectangle(box,tag='seqlabel',fill='yellow')
                    c.lift(obj)

        return

    def remove_seq_label(self, event):
		"""Remove sequence label"""
		c=self.seqframe
		c.delete('seqlabel')
		return

    def update_sequence_window(self,junk=None):
        """Delete everything from last time if needed"""
        for obj in self.seq_win_objs.keys():
            self.seqframe.delete(obj)
            del self.seq_win_objs[obj]

        # Do we have a DNA sequence?
        if not self.data['DNAseq']:
            return
          
        # Set base scaling and font        
        self.base_scale = self.base_scale_input.get()
        fontsize = self.seqfontsize_input.get()
        fontstyle = 'normal'
        if self.fontstyle_input.get() == 0:
            fontstyle = 'normal'
        elif self.fontstyle_input.get() == 1:
            fontstyle = 'bold'
        elif self.fontstyle_input.get() == 2:
            fontstyle = 'italic'

        self.restrfont = self.restr_site_font_input.get()+" 10"
        #create a font object this time - more convenient
        if fontstyle != 'italic':
            self.seqfont = tkFont.Font (family=self.seqfont_input.get(), size=fontsize,
                                    weight=fontstyle)
        else:
            self.seqfont = tkFont.Font (family=self.seqfont_input.get(), size=fontsize,
                                    slant="italic")

        # Change y scale to prevent crowding for large font sizes
        # this is done by getting the font size and applying a scale factor
        self.y_scale=1  #reset
        if fontsize>10:
            self.y_scale = self.seqfontsize_input.get()/10.0
            if self.y_scale>1:
                 self.y_scale=pow(self.y_scale,0.05)

        max_x,max_y=self.get_base_pos_on_screen(len(self.data['DNAseq']))
        max_x=max_x+30
        if max_x>self.canvas_x:
            #
            # Reconfigure
            #
            self.canvas_x=max_x
            self.canvas_y=self.canvas_height
            self.seqframe.configure(scrollregion=(0,0,self.canvas_x,self.canvas_y))
        #
        # Print the sequence in the window
        #
        self.seq_win_objs[self.seqframe.create_text(0,self.seq_row,font=self.labelfont,text='DNA seq',
            fill='darkgreen',anchor='w')]=1
        count=0
        colour='#999900'
        #print 'COLOUR= ',self.colour_seq_var.get()
        for letter in self.data['DNAseq']:
            if self.colour_seq_var.get()==0:
                colour='#66CC00'
            elif self.colour_seq_var.get()==1:
                cls={'A':'blue','T':'red','G':'darkgreen','C':'orange'}
                colour=cls[letter]
            elif self.colour_seq_var.get()==2:
                if count % 3 == 0 and colour == '#66CC00':
                    colour='#999900'
                elif count % 3 == 0 and colour=='#999900':
                    colour='#66CC00'
            else:
                colour='#66CC00'
            count=count+1
            x,y=self.get_base_pos_on_screen(count)
            self.seq_win_objs[self.seqframe.create_text(x,y,text=letter,font=self.seqfont,fill=colour,
                anchor='w')]=1

        # Print base sequence numbers
        position=0
        x,y=self.get_base_pos_on_screen(1)

        # Calculate y value once since it does not change - faster
        # If we have selected an ORF then we lift the base numbers a bit
        if self.data.has_key('ORF_selected'):
            y=(y-15)/self.y_scale

        # If a comparison sequence is displayed,
        # also raise the base number

        if self.show_comp_sequence.get()==1:
            if self.maxseqlevel==0:
                y=(y-15)/self.y_scale
            else:
                y=y-(self.maxseqlevel*15)/self.y_scale
            if self.primer_displayed==1 and self.maxseqlevel>0:
                y=(y-15)/self.y_scale
        else:
            y=(y-15)*self.y_scale

        for base in self.data['DNAseq']:
            position=position+1
            if float(position)/10.0==float(int(position/10)):
                x,ytmp=self.get_base_pos_on_screen(position+1)

                # Print them
                self.seq_win_objs[self.seqframe.create_line(x+3,y-5,x+3,y-20,fill='red',width=2)]=1
                self.seq_win_objs[self.seqframe.create_text(x+6,y-20,text='%d' %(position),font=self.labelfont,anchor='w')]=1

        # Translate in three forward frames
        import mutation
        AA_seqs3,AA_seqs1=mutation.translate(self.data['DNAseq'])
        self.data['AA_seqs1']=AA_seqs1
        self.data['AA_seqs3']=AA_seqs3
        #
        # Print the sequence
        #
        if self.data.has_key('ORF_selected'):
            self.orf_offset=8
            
            # Update the sequence in the ORF_selected array
            
            frame_sel=self.data['ORF_selected']['frame']
            start=self.data['ORF_selected']['start']
            stop=self.data['ORF_selected']['stop']
            self.data['ORF_selected']['aaseq']=self.data['AA_seqs1'][frame_sel-1][start-1:stop]
            self.data['ORF_selected']['aaseq3']=self.data['AA_seqs3'][frame_sel-1][start-1:stop]

            # Print only the coding sequence if we've selected an ORF
            seq=self.data['ORF_selected']['aaseq3']

            # Print the text
            x,y=self.get_aa_pos_on_screen(position,0)
            y=(y-self.orf_offset)*self.y_scale

            self.seq_win_objs[self.seqframe.create_text(0,y,text='ORF',font=self.labelfont,anchor='w')]=1

            # Backwards compatability
            if not self.data['ORF_selected']['aastart_number']:
                self.data['ORF_selected']['aastart_number']=1

            # Highlight the AA we're mutating
            highlight_aa=self.resnum.get()+self.data['ORF_selected']['start']-(self.data['ORF_selected']['aastart_number'])
            if highlight_aa<0:
                highlight_aa=None

            # Print sequence
            ORF_data=self.data['ORF_selected']
            self.print_aa3_seq(seq,y,
                               ORF_data['start']-1,ORF_data['frame'],
                               highlight_aa=highlight_aa)

            # Print the numbers
            self.print_aa_numbers(seq,15)
        else:
            # Print all three forward frames
            frame=0
            for AA_seq in AA_seqs1:
                x,y=self.get_aa_pos_on_screen(position,frame)
                y=y*self.y_scale
                self.seq_win_objs[self.seqframe.create_text(0,y,text='Frame %d' %(frame+1),
                font=self.labelfont,anchor='w')]=1
                self.print_aa_seq(AA_seq,y,0,frame)
                frame=frame+1
            self.print_aa_numbers(AA_seqs3[0])

        # Restriction digest
        # We delete the restriction details separately first, then re-plot
        self.clear_restriction_details()
        self.restriction_digest()

        # Overview window?
        self.update_overview_win()
        #move view to middle of scrollregion
        self.seqframe.yview("moveto", 0.2)
        self.seqframe.configure(bg=self.backgrcolor)
     
        return



    def print_aa_seq(self,sequence,y,position=0,frame=0):
        """
        Write a single aa sequence in the sequence window
        """
        for letter in sequence:
            position=position+1
            x,y_junk=self.get_aa_pos_on_screen(position,frame)
            self.seq_win_objs[self.seqframe.create_text(x,y,text=letter,font=self.seqfont,anchor='w')]=1
        return



    def print_aa3_seq(self,sequence,y,position=0,frame=0,highlight_aa=None):
        """
        Write a single aa sequence in the sequence window
        """
        for aa_name in sequence:
            position=position+1
            x,y_junk=self.get_aa_pos_on_screen(position-1,frame+1)
            if position==highlight_aa:
                colour='red'
            else:
                colour='black'
            self.seq_win_objs[self.seqframe.create_text(x+2,y,text=aa_name,font=self.seqfont,anchor='w',fill=colour)]=1
        return


    def print_aa_numbers(self,aa_seq,y_offset=55,frame=1):
        """
        Print amino acid sequence numbers
        """
        # First the label
        x,y=self.get_aa_pos_on_screen(0,frame)
        #determine y position based on y scaling - font size dependent
        y=(y+y_offset-self.orf_offset)*frame*self.y_scale

        self.seq_win_objs[self.seqframe.create_text(0,y,text='AA number',font=self.labelfont,anchor='w')]=1

        # Now the numbers
        position=0
        number_start=0
        offset=2
        if self.data.has_key('ORF_selected'):
            number_start=self.data['ORF_selected']['aastart_number']
        for aa in aa_seq:
            position=position+1
            #we skip zero by resetting the offset
            if position + number_start-offset == 0:
                offset=1
            aa_number=number_start+position-offset
            #print 'aa_number', aa_number, 'pos', position
            if aa_number == 1.0:
                self.seq_win_objs[self.seqframe.create_line(x+4,y,x+4,y-15,fill='red',width=2)]=1
            if float(aa_number)/5.0 == float(int(aa_number/5)) and aa_number!=0.0:
                if self.data.has_key('ORF_selected'):
                    ORF_data=self.data['ORF_selected']
                    x,y=self.get_aa_pos_on_screen(position+ORF_data['start']-1,frame)
                else:
                    x,y=self.get_aa_pos_on_screen(position,frame)
                y=(y+y_offset-self.orf_offset)*frame*self.y_scale
                self.seq_win_objs[self.seqframe.create_line(x+4,y,x+4,y-15,fill='red',width=2)]=1
                self.seq_win_objs[self.seqframe.create_text(x+7,y,text='%d' %(aa_number),font=self.labelfont,anchor='w')]=1
        return


    def get_base_pos_on_screen(self,position):
        """return the x and y position of the base on the screen"""

        #return self.seq_xstart+float(position-1)*8,self.seq_row
        return self.seq_xstart+float(position-1)*self.base_scale,self.seq_row

    def get_DNApos_fromcoords(self,x,y):
        """From X and Y coordinate, return the DNA base number"""

        # Are we close to the DNA sequence?
        if abs(y-self.seq_row)>10:
            return None

        # ok, DNA it is
        pos=int(float(x-self.seq_xstart+4.0)/self.base_scale)
        return pos


    def get_aa_pos_on_screen(self,position,frame):
        """
        # return the x and y position of the aa on the screen
        # frame is 0, 1 or 2
        """
        position=position*3+float(frame)-1
        x,y=self.get_base_pos_on_screen(position)
        y=y+20.0+float(frame)*15.0
        return x,y

    #
    # Menu actions
    #

    def invert_seq(self):
        """
        Invert the DNA sequence
        """
        if not self.data['DNAseq']:
            self.invert_seq_var.set(0)
            self.warning('No DNA sequence loaded','You have to load a DNA sequence first')
            return
        inverted=''
        for count in range(len(self.data['DNAseq'])):
            pos=-count-1
            inverted=inverted+self.data['DNAseq'][pos]
        self.data['DNAseq']=inverted
        #
        # Update
        #
        self.update_sequence_window()
        return


    def complementary_seq(self):
        """
        Get the complementary DNA sequence
        """
        if not self.data['DNAseq']:
            self.complement_seq_var.set(0)
            self.warning('No DNA sequence loaded','You have to load a DNA sequence first')
            return
        compl={'A':'T','T':'A','C':'G','G':'C'}
        comDNA=''
        for base in self.data['DNAseq']:
            comDNA=comDNA+compl[base]
        self.data['DNAseq']=comDNA
        #
        # Update
        #
        self.update_sequence_window()
        return


    def warning(self,text1,text2):
        """
        Show a warning
        """
        import tkMessageBox
        tkMessageBox.showwarning(text1,text2)
        return


    def quit(self,event=None):
        """
        Exit or return to DB_Main.
        If going to DB_Main, make sure that we save the primer_dict in the Database.
        Make sure that we close the primer database window.
        """

        #Make sure that we close the primer database window if open
        if self.pDB_open:
            self.pDB_open.close_pDB()
            self.pDB_open=None
        if self.sequ_win:
            self.sequ_win.close()
            self.sequ_win=None

        # Close the design_mutagenic primer window if open
        if getattr(self,'primer_win',None):
            self.primer_win.destroy()

        # Manage the data transfer - or simply exit
        if self.parent:
            # Check if the altered flag is set
            if self.data['DNAseq_status']=='ALTERED':
                # DNA sequence is altered - we have to ask if a new PEAT_DB record should be created

                import tkMessageBox
                if not tkMessageBox.askyesno('DNA sequence altered',
                                             'The DNA sequence has been altered.'
                                             'Do you want to create a new mutant record?',
                                             parent=self.master):

                    # Show another warning
                    self.data['DNAseq_status']='MODIFY RECORD'
                else:

                    # Create a new PEAT_DB record
                    self.data['DNAseq_status']='ADD NEW RECORD'

            elif self.data['DNAseq_status']!='ALTERED' and self.data['DNAseq_status']=='ORF_CHANGED' and self.data.has_key('ORF_selected'):
                # Only the ORF is altered, ask if we want to save
                import tkMessageBox
                if tkMessageBox.askyesno('ORF shifted',
                                         'The ORF has been moved. Do you want to try to modify the record?',
                                          parent=self.master):
                    self.data['DNAseq_status']='MODIFY ORF'

            elif self.data['DNAseq_status']=='NEW SEQ':
                #this only happens when it has no previous sequence
                pass

            import copy
            if self.parent.DB.data:
            #if self.parent.data['DBinstance']:

                # Record the changes in the primer database
                # First delete everything in the old parent version of the primer db

                self.parent.DB.meta['DNAtool_primers']={}
                # Now add the current primers to the parent instance of the db
                for primer in self.data['primer_dict'].keys():
                    self.parent.DB.meta['DNAtool_primers'][primer]={'sequence':self.data['primer_dict'][primer]['sequence']}
                    self.parent.DB.meta['DNAtool_primers'][primer]['description']=self.data['primer_dict'][primer]['description']

                    if not self.data['primer_dict'][primer].has_key('description'):
                        self.data['primer_dict'][primer]['description']='None entered'

                del self.data['primer_dict']
                self.parent.DNAtool_state = copy.deepcopy(self.data)

            self.master.destroy()
        else:
            import os
            os._exit(0)
        return



    def seq_display_settings(self):
        """Dialog to change sequence display settings"""
        # Open a new window for setting the restriction enzymes

        self.seq_display_setupwin=Toplevel()
        self.seq_display_setupwin.geometry('+300+450')
        self.seq_display_setupwin.title('Sequence Display Setup')

        # Spacing between bases
        row=1
        lblspace=Label(self.seq_display_setupwin,text='Bases Spacing:')
        lblspace.grid(row=row,column=0,padx=3,pady=2)
        bscaleentry=Scale(self.seq_display_setupwin,from_=8,to=20,resolution=1,orient='horizontal',
                            relief='ridge',variable=self.base_scale_input,label='scale factor')
        bscaleentry.grid(row=row,column=1, sticky='wens', padx=3,pady=2)
        row=2
        lblfont=Label(self.seq_display_setupwin,text='Seq Font:')
        lblfont.grid(row=row,column=0,padx=3,pady=2)
        fontentry_button=Menubutton(self.seq_display_setupwin,textvariable=self.seqfont_input,
					relief=RAISED,width=16)
        restr_fontentry_button=Menubutton(self.seq_display_setupwin,textvariable=self.restr_site_font_input,
					relief=RAISED,width=16)
        fontentry_menu=Menu(fontentry_button,tearoff=0)
        restr_fontentry_menu=Menu(restr_fontentry_button,tearoff=0)
        fontentry_button['menu']=fontentry_menu
        restr_fontentry_button['menu']=restr_fontentry_menu

        # Other fonts available
        fts=['Arial','Courier','Verdana','Fixed','Times']
        for text in fts:
            #text='Font '+text
            fontentry_menu.add_radiobutton(label=text,
                                            variable=self.seqfont_input,
                                            value=text,
                                            indicatoron=1)
            restr_fontentry_menu.add_radiobutton(label=text,
                                            variable=self.restr_site_font_input,
                                            value=text,
                                            indicatoron=1)
        fontentry_button.grid(row=row,column=1, sticky='nes', padx=3,pady=2)

        row=3
        lblfontsize=Label(self.seq_display_setupwin,text='Sequence Font Size:')
        lblfontsize.grid(row=row,column=0,padx=3,pady=2)
        fontsizeentry=Scale(self.seq_display_setupwin,from_=8,to=20,resolution=1,orient='horizontal',
                            relief='ridge',variable=self.seqfontsize_input)

        fontsizeentry.grid(row=row,column=1, sticky='wens',padx=3,pady=2)
        row=4
        frame = Frame(self.seq_display_setupwin)
        fontstyle_label = Label(frame, text='Font Style:')
        fontstyle_label.grid(row=0,column=0)
        fontstyle = Radiobutton(frame, text="plain", variable=self.fontstyle_input, value=0)
        fontstyle1 = Radiobutton(frame, text="bold", variable=self.fontstyle_input, value=1)
        fontstyle2 = Radiobutton(frame, text="italic", variable=self.fontstyle_input, value=2)
        fontstyle.grid(row=0,column=1)
        fontstyle1.grid(row=0,column=2)
        fontstyle2.grid(row=0,column=3)
        frame.grid(row=row,column=0,columnspan=2,sticky='news', padx=3,pady=2)

        row=5
        self.backgrcolorbutton = Button(self.seq_display_setupwin, text='background color', bg=self.backgrcolor,
                                    command=self.setbackgrcolor)
        self.backgrcolorbutton.grid(row=row,column=1, sticky='nes', padx=3,pady=2)
        row=6
        restrfont=Label(self.seq_display_setupwin,text='Restr. Site Font:')
        restrfont.grid(row=row,column=0,padx=3,pady=2)
        restr_fontentry_button.grid(row=row,column=1, sticky='nes', padx=3,pady=2)

        row=7
        
        # Apply Button        
        b = Button(self.seq_display_setupwin, text="Apply Settings", command=self.update_window_formatting)
        b.grid(row=row,column=1,sticky='wens',padx=4,pady=4)
        
        # Close button    
        c=Button(self.seq_display_setupwin,text='Close',command=self.close_seq_display_setupwin)
        c.grid(row=row,column=0,sticky='wens',padx=4,pady=4)
        
        # Save Settings button        
        row=8
        c=Button(self.seq_display_setupwin,text='Save as Default',command=self.save_preferences)
        c.grid(row=row,column=0,columnspan=2,sticky='wens',padx=4,pady=4)
        return

    def setbackgrcolor(self):
        clr = self.getaColor(self.backgrcolor)
        if clr != None:
            self.backgrcolor = clr
            self.backgrcolorbutton.configure(bg=clr)
        return
    
    def getaColor(self, oldcolor):
        import tkColorChooser
        ctuple, newcolor = tkColorChooser.askcolor(title='pick a color', initialcolor=oldcolor,
                                                   parent=self.seq_display_setupwin)
        if ctuple == None:
            return None
        return str(newcolor)

    def close_seq_display_setupwin(self):
        """Closes the display settings window"""
        self.seq_display_setupwin.destroy()
        return


    def save_preferences(self):
        """Saves the sequence display settings to prefs file"""

        print 'Saving DNAtool preferences'
        self.preferences.set('seqfont',self.seqfont_input.get())
        self.preferences.set('seqfontsize',self.seqfontsize_input.get())
        self.preferences.set('fontstyle',self.fontstyle_input.get())
        self.preferences.set('base_scale',self.base_scale_input.get())
        self.preferences.set('restr_font',self.restr_site_font_input.get())
        self.preferences.set('backgrcolor',self.backgrcolor)
        print self.preferences.get('restr_font')
        return

    def load_preferences(self):
        """Loads the sequence display settings from prefs file, if present"""

        print 'Loading current DNAtool preferences'
        self.preferences=Preferences('DNAtool',{'canvas_height':600})

        try:
            f=self.preferences.get('seqfont')
            self.seqfont_input.set(f)
        except:
            self.preferences.set('seqfont','Courier')
        try:
            f=self.preferences.get('seqfontsize')
            self.seqfontsize_input.set(f)
        except:
            self.preferences.set('seqfontsize',12)
        try:
            f=self.preferences.get('fontstyle')
            self.fontstyle_input.set(f)
        except:
            self.preferences.set('fontstyle',1)
        try:
            f=self.preferences.get('base_scale')
            self.base_scale_input.set(f)
        except:
            self.preferences.set('base_scale',10)
        try:
            f=self.preferences.get('restr_font')
            self.restr_site_font_input.set(f)
        except:
            self.preferences.set('restr_font','Arial')
        try:
            self.backgrcolor = self.preferences.get('backgrcolor')
        except:
            pass                
        return

            
    def update_window_formatting(self):
        """Do this only when display setup changes are applied. Ensures all objects, including
           applied primers are updated with the new formatting."""
        self.update_sequence_window()
        if self.pDB_open:
           self.pDB_open.refresh_primer()
        if self.show_comp_sequence.get==1:
           self.sequ_win.refresh_DNAseq()
        return


    def overview_on_off(self):        
        """Turn the over view window on and off"""
        
        if self.overview_win:
            self.overview_button.deselect()
            self.overview_win.destroy()
            self.overview_win=None
        else:
            self.overview_button.select()
            if not self.data.has_key('AA_seqs1'):
                self.warning('No DNA sequence loaded','Load a DNA sequence first')
                self.overview_button.deselect()
                return
            #
            # Open Canvas and draw lines
            #
            self.overview_win=Toplevel()
            self.overview_win.geometry('300x100+400+350')
            self.overview_win.title('Open reading frames')
            self.overview_frame=Canvas(self.overview_win,bd=5,bg='white',width=300,height=150)
            self.overview_frame.xview("moveto", 0)
            self.overview_frame.yview("moveto", 0.2)
            self.overview_frame.grid(row=0,column=0)
            #
            # Draw
            #
            self.update_overview_win()
        return


    def open_pDB(self, overwrite=False):
        """Open the primer database window"""
        import primer_database
        print self.pDB_open
        if not self.pDB_open or overwrite == True:
            if self.pDB_open:
                self.pDB_open.close_pDB()
            self.pDB_open=primer_database.primer_database(self, parentframe=self.master_win)
        return

    def open_primer_design(self,event=None):
        if not self.pDB_open:
            self.design_mutagenic_primer(self)
            print 'calling primer design without pdb window open'
        else:
            self.design_mutagenic_primer(self.pDB_open)
            print 'calling primer design with pdb window open'
        return

    def update_overview_win(self):
        """
        Update overview
        """
        # Is the window open?
        if not self.overview_win:
            return

        # Delete old objects
        for obj in self.overview_win_objs.keys():
            self.overview_frame.delete(obj)
            del self.overview_win_objs[obj]

        # Update the overview window
        frame=0
        for AA_seq in self.data['AA_seqs1']:
            position=0
            x=60
            y=frame*25+25
            self.overview_win_objs[self.overview_frame.create_text(0,y,text='Frame %d' %(frame+1),font=self.seqfont,anchor='w')]=1
            self.overview_win_objs[self.overview_frame.create_line(x,y,x+200,y,fill='green',width=2)]=1
            for letter in AA_seq:
                position=position+1
                xpos=x+(200.0/float(len(AA_seq))*float(position))
                #x,y=self.get_aa_pos_on_screen(position,frame)
                if letter=='*':
                    self.overview_win_objs[self.overview_frame.create_line(xpos,y,xpos,y-10,fill='red',width=2)]=1
                    self.overview_win_objs[self.overview_frame.create_text(xpos,y-15,text='%d' %(position),font='Times 6',anchor='w')]=1
            frame=frame+1
        return

    def select_ORF(self):
        """
        Clear the ORF_selected flag
        """
        if self.data.has_key('ORF_selected'):
            del self.data['ORF_selected']
            self.data['DNAseq_status']='ORF_CHANGED'
        #if no orf previously then this must be a new record?
        else:
            self.data['DNAseq_status']='NEW SEQ'

        # Start the ORF selection
        import ORF_module
        X=ORF_module.ORF_handler(self)
        if self.data.has_key('ORF_selected'):
            start=self.data['ORF_selected']['start']
            stop=self.data['ORF_selected']['stop']
            frame=self.data['ORF_selected']['frame']
            ORF_num=self.data['ORF_selected']['number']
            #
            base_start=(start-1)*3+(frame-1)
            base_stop=stop*3+(frame-1)

            # mark that we have an ORF and get rid of all useless info
            self.data['Project saved']=None
            #
            # Activate the primer design button(s)
            #
        else:
            pass
        #
        # Update everything
        #
        self.assess_status()
        # Update view
        self.update_sequence_window()
        #self.data['DNAseq_status']='ALTERED'

        return


    def about_DNAtool(self):
        """Display about window"""
        self.ab_win=Toplevel()
        self.ab_win.geometry('+100+350')
        self.ab_win.title('About DNAtool')
        row=0
        self.logo=Text(self.ab_win)
        imgpath = os.getcwd()
        imgfile = os.path.join(imgpath,'DNAtool.gif')
        print imgfile

        photo=PhotoImage(file=imgfile)
        label = Label(self.ab_win,image=photo)
        label.image = photo # keep a reference!
        #label.pack()
        label.grid(row=row,column=0,sticky='news')
        #import base64
        #print "icon='''\\\n" + base64.encodestring(open(imgfile, "rb").read()) + "'''"

        text=['DNAtool ','Component of Protein Engineering and Analysis Tool',
        'A stand-alone subcomponent of PEAT that contains all standard functions',
        'for DNA sequence analysis, designed for the researcher who works',
        'with smallish pieces of DNA containing a single gene.',
        'Authors: Jens Erik Nielsen & Damien Farrell, University College Dublin',
        '(C) Copyright 2003- Jens Erik Nielsen All rights reserved']

        for line in text:
            tmp=Label(self.ab_win,text=line)
            tmp.grid(row=row,column=0,sticky='news')
            row=row+1
        return


    def online_documentation(self,event=None):
        """Open the online documentation"""
        import webbrowser
        link='http://enzyme.ucd.ie/PEAT/'
        webbrowser.open(link,autoraise=1)
        return


    def compare_DNA_seq(self):
        """load and display a sequence from a file to compare to the currently
            displayed sequence"""

        import tkMessageBox
        if self.data['DNAseq']=='' or not self.data['DNAseq']:
            tkMessageBox.showwarning(
            "Compare DNA Sequence",
            "Load a DNA sequence first\nUse the Browse Button"
            )
            return
        else:
            self.show_comp_sequence.set(1)
            import DNA_sequencing
            if not self.sequ_win:
                self.sequ_win=DNA_sequencing.sequencing_window(self)
            self.sequ_win.show_DNAseq()
            if self.pDB_open:
                self.pDB_open.refresh_primer()

        return

    #slight duplication with this function, would be better to call both
    #from one function with different arguments
    def compare_multiple_seqs(self):
        """load and display many sequences in a set of files or zip file"""
        import tkMessageBox
        if self.data['DNAseq']=='' or not self.data['DNAseq']:
            tkMessageBox.showwarning(
            "Compare DNA Sequence",
            "Load a DNA sequence first\nUse the Browse Button"
            )
            return
        else:
            self.show_comp_sequence.set(1)
            import DNA_sequencing
            if not self.sequ_win:
                #print self,'inDNAtool'
                self.sequ_win=DNA_sequencing.sequencing_window(self)
            self.sequ_win.show_multiple_DNAseq_dialog()
            #self.wait_window(self.sequ_win.si_window)
        return


    def clear_comp_seqs(self, event=None):
        """switch off display of extra comparison sequences"""

        self.show_comp_sequence.set(0)
        print 'clearing comp seq'
        self.seqframe.delete('comparison_seq')
        self.seqframe.delete('strip')
        self.max_seq_level=0
        self.update_sequence_window()
        if self.pDB_open:
            self.pDB_open.refresh_primer()
        return


if __name__=="__main__":
    """Run the app"""
    import sys, os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="projectfile",
                        help="Open a project", metavar="FILE")

    opts, remainder = parser.parse_args()

    if opts.projectfile != None:
        import pickle
        fd=open(opts.projectfile)
        data=pickle.load(fd)
        fd.close()
        app = MainWindow(data_passed=data)
    else:
        app = MainWindow()

    app.mainloop()



