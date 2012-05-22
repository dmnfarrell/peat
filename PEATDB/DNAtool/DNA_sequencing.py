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
import tkFileDialog
import Pmw
import os, zipfile, types
import evaluate_primer
import primer_alignment
import PEATDB.ProgressBar as ProgressBar

class sequencing_window(primer_alignment.align_primer):

    def __init__(self,main_self):

        self.main_self = main_self
        self.parent = main_self        
        self.choose_seq_win=None
        self.main_self.maxseqlevel=0
        return

    #
    # Base functions for load and display of a sequence
    # This set of classes has similar functionality to the primer db gui
    # The gui dialog is provided in a seperate function below
    #
    def show_DNAseq(self,seqfile=None,usecurrent=0,clearprev=1,level=0):
        """Load and display a sequence from a file"""

        if seqfile==None:
            seqfile = self.main_self.dnaseq_read(newprotein=None,variable='seq_comp',var2='seq_comp_file')
        # Pass the file to dnaseq_read instead of getting it to open a dialog
        else:
            #print 'seqfile:',seqfile
            self.main_self.dnaseq_read(newprotein=None,fileinput=seqfile,variable='seq_comp',
                                        var2='seq_comp_file')
        #print self.main_self.data['seq_comp']
        # Clear previous seq first, this is default value in parameter list
        if clearprev==1:
            self.main_self.seqframe.delete('comparison_seq')
        #print self.main_self.data['seq_comp']
        #print
        #print self.main_self.data['DNAseq']

        #if usecurrent is 1, skip the calculation and use current seq (only for refresh)
        if usecurrent==0:
            sites=self.find_primer_binding_sites(self.main_self.data['seq_comp'],
                                                    self.main_self.data['DNAseq'])
            print sites
            #
            # Print the number of binding sites
            #
            best_score=0
            best_position=None
            scores=[]
            first_neg=1
            for position,score in sites:
                scores.append(score)
                #
                # Find the best position
                #
                if score>best_score:
                    best_score=score
                    best_position=position
            this_seq={'startpos':best_position}
            this_seq['sequence']=self.main_self.data['seq_comp']
            self.current_seq=this_seq

        if self.main_self.show_comp_sequence.get()==1:            
            x = os.path.split(seqfile)
            self.draw_primer(self.current_seq,colour='#800080',thetag='comparison_seq',seqname=x[1],lvl=level)
        return


    def show_multiple_DNAseq_dialog(self):
        """
        This function provides a dialog to view and select the multiple loaded sequences
        """
        self.choose_seq_win=Toplevel()
        self.choose_seq_win.title("Load Multiple DNA Sequences")
        #self.choose_seq_win.geometry("%dx%d%+d%+d" % (40, 20, 0, 0))
        self.choose_seq_win.geometry("+%d+%d" % (self.parent.master_win.winfo_rootx()+550,
                                  self.parent.master_win.winfo_rooty()+100))

        OpenSeqFilesbtn=Button(self.choose_seq_win,text='Load Files',command=self.load_multiple)
        OpenSeqFilesbtn.grid(row=0,column=0,sticky='news',padx=3,pady=3)
        OpenSeqZipbtn=Button(self.choose_seq_win,text='Load Zip',command=self.load_Zip)
        OpenSeqZipbtn.grid(row=0,column=1,sticky='news',padx=3,pady=3)
        Removeitemsbtn=Button(self.choose_seq_win,text='Remove All',command=self.clear_list)
        Removeitemsbtn.grid(row=0,column=2,sticky='news',padx=3,pady=3)

        import tkFont
        detailsfont = tkFont.Font ( family="Helvetica", size=10, weight="bold" )
        yscrollbar=Scrollbar(self.choose_seq_win,orient='vertical',width=11)
        yscrollbar.grid(row=1,column=4,sticky='news',padx=2)

        self.details=Listbox(self.choose_seq_win,bg='white',
                             fg='black',
                             height=15,
                             width=40,
                             yscrollcommand=yscrollbar.set,
                             font=detailsfont,
                             selectmode=EXTENDED)

        yscrollbar.config(command=self.details.yview)
        self.details.grid(row=1,column=0,columnspan=3,sticky='NWS',padx=3,pady=3)
        self.details.config(state=NORMAL)

        Clear_CloseSeqbtn=Button(self.choose_seq_win,text='Close',command=self.close)
        Clear_CloseSeqbtn.grid(row=3,column=0,sticky='news',padx=4,pady=3)
        Clearbtn=Button(self.choose_seq_win,text='Clear Display',command=self.main_self.clear_comp_seqs)
        Clearbtn.grid(row=3,column=1,columnspan=2,sticky='news',padx=4,pady=3)

        self.seq_progressbar = ProgressBar.ProgressBar(self.choose_seq_win)
        self.seq_progressbar.frame.grid(row=4,column=1,columnspan=2,padx=2,pady=4)
        progrlbl = Label(self.choose_seq_win,text='Progress:')
        progrlbl.grid(row=4,column=0,sticky='news',padx=2,pady=4)
        self.choose_seq_win.grid_rowconfigure(1, weight=1)

        #do balloon help for buttons
        try:
            import Pmw
        except:
            import tkMessageBox
            tkMessageBox.showinfo("Missing component","To run this application you need to install Pmw (Python Megawidgets)\n or do 'Force update' to get the latest version of PEAT which includes Pmw'.")
        self.balloon=Pmw.Balloon(self.choose_seq_win)
        help=[[OpenSeqFilesbtn,'Load a sequence from a list of files'],
              [OpenSeqZipbtn,'Open from a zip file'],
              [Removeitemsbtn,'Clear list of sequences'],
              [Clear_CloseSeqbtn,'Clear sequences from the display and close dialog'],
              [Clearbtn,'Clear sequences from the display']]
        for btn,txt in help:
            self.balloon.bind(btn,txt)

        #self.details.bind("<Button-1>",self.display_current)
        self.details.bind("<ButtonRelease-1>",self.display_current)
        self.details.bind("<Double-Button-1>",self.change_seq_direction)
        return


    def refresh_DNAseq(self):
        """Redisplay the currently shown sequence in the main window"""
        self.show_DNAseq(usecurrent=1)
        return

    # Load multiple files and show them in the details listbox, put the first one
    # on display in the sequence window
    #
    def load_multiple(self,event=None):
        """load multiple files and show them in the details widget"""
        # a list with the current sequences loaded for display
        self.seq_list=[]
        self.seq_list_direction={}
        self.clear_list()
        # populate file list here
        self.seq_list = self.open_filelist()
        for name in self.seq_list:           
            x = os.path.split(name)
            self.details.insert(END,x[1])
            self.seq_list_direction[name]='f'

        return

    # Display the current sequence(s) selected in the listbox
    def display_current(self,event=None):
        """Display the current sequences selected in the listbox"""
       
        # Reset this if main window has been cleared while listbox open
        self.main_self.show_comp_sequence.set(1)
        self.main_self.maxseqlevel=0
        self.main_self.seqframe.delete('strip')
        # Figure out the selected seq and display in textbox
        #tmp_sel=self.details.curselection()
        tmp_sel = map(int, self.details.curselection())
        p=0 # set positions above main seq
        if len(tmp_sel)==1:
            num=tmp_sel[0]
            self.show_DNAseq(self.seq_list[num])
        else:
            # Display multiple sequences at once by creating a current list
            self.main_self.seqframe.delete('comparison_seq')
            self.current_list=[]

            max=len(tmp_sel)
            #print 'max',max
            for num in tmp_sel:
                self.current_list.append(self.seq_list[num])
                self.show_DNAseq(self.seq_list[num],clearprev=0,level=p)
                p=p+1
                #store the current highest sequence position in the parent instance
                #this is so any primer can be displayed above that
                self.main_self.maxseqlevel=p
                #print 'P=',p
                m=(p/float(max))*100
                #print 'amount done',m
                self.seq_progressbar.updateProgress(newValue=m)

            #print self.current_list
        if p>0:
            self.draw_seq_boxes(p)
        #refresh the primer if present
        if self.main_self.pDB_open:
            self.main_self.pDB_open.refresh_primer()


    def change_seq_direction(self,event=None):
        """Sets forward or reverse sense for entry based on single/double-click"""
        #pop up dialog box here first time this entry is selected
        #sel=self.details.curselection()
        items = map(int, self.details.curselection())
        num=items[0]
        name=self.seq_list[num]
        if self.seq_list_direction[name]=='f':
            self.details.itemconfig(num,fg='red')
            self.seq_list_direction[name]=='r'
        else:
            self.details.itemconfig(num,fg='black')
            self.seq_list_direction[name]=='f'
        return

    def load_Zip(self,event=None):
        """Loads sequences from a zip file and show items in the details widget"""

        zfilename = None
        zfilename = self.open_zipfile()
        if zfilename == None:
            return
        zfile = zipfile.ZipFile( zfilename, "r" )
        zfile.printdir()
        self.seq_list=[]
        for info in zfile.infolist():
            fname = info.filename
            if fname.endswith(".clipped"):
                self.seq_list.append(fname)
            # decompress each file's data
            #data = zfile.read(fname)
            #print data
                x = os.path.split(fname)
                self.details.insert(END,x[1])
        return

    # Clears the current sequences from the listbox
    def clear_list(self,event=None):
        """Clear listbox"""
        self.details.delete(0, END)
        return

    #
    # Draw a grey shaded boxes to highlight all displayed sequences
    #

    def draw_seq_boxes(self,levels=0):

        c=self.main_self.seqframe
        #
        # Delete all the old boxes first
        bgcolour='gray10'
        c.delete('strip')

        # Draw the box
        x,y=self.parent.get_base_pos_on_screen(0)
        x=0
        y=y-7.5
        height=15
        end=self.main_self.canvas_x
        for l in range(levels):
            print 'l:',l%2
            if l%2==0:
                bgcolour='gray95'
            else:
                bgcolour='gray90'
            obj=c.create_rectangle( x,y, x+end,y-height,
                                    fill=bgcolour,width=0,
                                    tag='strip')
            y=y-height
            c.tag_lower(obj)
        return


    def close(self,event=None):
        """Closes the multiple sequences window"""

        self.main_self.clear_comp_seqs()
        if self.choose_seq_win:
            self.choose_seq_win.destroy()
        return

    # Open multiple files dialog
    def open_filelist(self):
        """Open multiple filenames list dialog"""
        
        filelist=tkFileDialog.askopenfilenames(defaultextension='.seq.clipped',
                                                initialdir=self.main_self.data['datadir'],
                                                filetypes=[("Clipped Seq","*.seq.clipped"),
                                                           ("PIR file","*.pir"),
                                                           ("FASTA file","*.txt"),
                                                           ("GenBank file","*.gb"),
                                                           ("BSML file",".xml"),
                                                           ("All files","*.*")],
                                                parent=self.parent.master)
        if type(filelist) == types.UnicodeType:
            #master_win is reference to parent root
            filelist = self.parent.master_win.tk.splitlist(filelist)
            
        return filelist

    # Open zip files dialog
    def open_zipfile(self):
    
        filename=tkFileDialog.askopenfilename(defaultextension='.zip',
                                                initialdir=self.main_self.data['datadir'],
                                                filetypes=[("Zip file","*.zip"),
                                                           ("All files","*.*")],
                                                parent=self.parent.master)
        return filename


