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

from Tkinter import *
import tkMessageBox, tkFileDialog
import Pmw

class multiple_calc_analysis:
    """Class for analysing multiple pKa calculations simultaneously"""

    def __init__(self,parent):
        """open the window and display all"""
        self.aas={'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','TPO':'Z','SEP':'X'}
     
        self.parent=parent
        self.multwin=Toplevel()
        self.parent.set_position(self.multwin)
        self.multwin.title('Multiple pKa calculation analysis')
        #
        # Canvas
        #
        Pmw.initialise()
        self.canv = Pmw.ScrolledCanvas(self.multwin,
                borderframe = 1,
                labelpos = 'n',
                label_text = 'pKa Sequence Display',
                usehullsize = 1,
                hull_width = 800,
                hull_height = 400,
                hscrollmode='dynamic',
                vscrollmode='dynamic'
        )
        self.canv.interior().configure(bg='white')
        self.canv.grid(row=0,column=0)

        self.multwin.update_idletasks()
        
        import copy
        self.multcalcs=copy.deepcopy(self.parent.calcs)
        self.get_sequences()

        self.colour_function=self.colour_dpKa
        self.display_sequences()
        self.canv.resizescrollregion()
        """Pulldown menu"""
        self.menu=Menu(self.multwin)     
        self.anal_menu=Menu(self.menu,tearoff=0)
        self.anal_menu.add_command(label='Set colour scheme',command=self.set_colours)
        self.anal_menu.add_command(label='Save all sequences',command=self.save_sequences)
        self.anal_menu.add_command(label='Import clustal aln file',command=self.get_clustal_aln)
        self.anal_menu.add_command(label='Align similar pKa shifts',command=self.align_pkas)
        self.anal_menu.add_command(label='Find pka signatures',command=self.find_signatures)
        self.anal_menu.add_command(label='Quit',command=self.multwin.destroy)
        self.menu.add_cascade(label='Analyse',menu=self.anal_menu)
        self.multwin.config(menu=self.menu)        
        return

    #
    # -----
    #

    def set_colours(self):
        self.colour_function=self.colour_dpKa
        return


    #
    # ------
    #


    def save_sequences(self):
        """Save all sequences in a single PIR file for an alignment program"""
        #
        # Just do it
        #
        import tkFileDialog, os, pickle
        savedir=os.getcwd()
        filename=tkFileDialog.asksaveasfilename(defaultextension=".pir",
                                                initialfile='sequences.pir',
                                                initialdir=savedir,
                                                filetypes=[("Pir file","*.pir")],
                                                parent=self.multwin)
        if filename:
            fd=open(filename,'w')
            for calc1 in self.multcalcs.keys():
                seq1=self.multcalcs[calc1]['aa1seq']
                fd.write('>P1;%s\n\n%s*\n\n' %(calc1,seq1))
            fd.close()
        return

    #
    # ----
    #

    def get_clustal_aln(self):
        """Load a clustal aln file and adjust the alignment of the sequences accordingly"""
        import tkFileDialog, os
        alnfilename=tkFileDialog.askopenfilename(defaultextension='.aln',
                                                 initialdir=os.getcwd(),
                                                 filetypes=[("Clustal alignment file","*.aln"),
                                                            ("All files","*.*")],
                                                 parent=self.multwin)
        if alnfilename:
            import EAT_DB.sequence_alignment as SA
            self.alignment=SA.read_clustal_aln_file(alnfilename)
            self.display_sequences()
        return

    #
    # ----
    #
    
    def colour_dpKa(self,calc,resnum,restype):
        """Colour as a function of dpKa"""
        resid=resnum+':'+restype
        if calc['pkavals'].has_key(resid):
            dpka=calc['pkavals'][resid]['pKa']
            if dpka:
                dpka=dpka-calc['pkavals'][resid]['modelpK']
            else:
                dpka=None
            return dpka 
        return 'nontit'


    def align_pkas(self):
        return

    def find_signatures(self):
        return

    def quit_multanal(self):
        return

    #
    # -----
    #

    def get_sequences(self):
        """Extract all sequences from the PDB files"""
        import Protool
        for calc in self.multcalcs.keys():
            X=Protool.structureIO()
            X.parsepdb(self.multcalcs[calc]['pdblines'])
            self.multcalcs[calc]['sequence']=X.sequence[:]
            aa1=''
            for number,restype in self.multcalcs[calc]['sequence']:
                aa1=aa1+self.aas[restype]
            self.multcalcs[calc]['aa1seq']=aa1
        return

    #
    # -----
    #

    def display_sequences(self):
        """Display all the sequences"""
        if not getattr(self,'sequencesobjs',None):
            self.sequencesobjs=[]
        #
        # Delete all graphics objects from last round
        #
        for obj in self.sequencesobjs:
            self.canv.delete(obj)
        self.sequencesobjs=[]
        #
        # Display the new sequences
        #
        calcs=self.multcalcs.keys()
        calcs.sort()
        seq_num=0
        for calc in calcs:
            seq_num=seq_num+1
            if getattr(self,'alignment',None):
                alignment=self.alignment[calc]
            else:
                alignment=None
            self.display_sequence(self.multcalcs[calc],calc,seq_num,alignment)
        return

    #
    # ----
    #

    def display_sequence(self,calculation,calculation_name,seq_num,alignedseq=None):
        """Display a single sequence in the window"""
        y_spacing=15
        y=20+seq_num*y_spacing
        x_start=100
        spacing=9
        count=0
        display_count=0
        aligncount=0
        self.canv.create_text(0,y,text=calculation_name,anchor='nw')
        for aanum,restype in calculation['sequence']:
            aa1=self.aas[restype]
            #
            # If we have an alignment then skip forward until we have a residue
            #
            if alignedseq:
                #print alignedseq[aligncount],aa1,count,aligncount
                if alignedseq[aligncount]==aa1:
                    pass
                else:
                    if alignedseq[aligncount]!='-':
                        #print alignedseq[aligncount]
                        #print alignedseq
                        print 'error in alignment'
                        #print calculation['sequence']
                        stop
                    else:
                        while alignedseq[aligncount]=='-':
                            x=x_start+display_count*spacing
                            self.sequencesobjs.append(self.canv.create_text(x,y,text='-',fill='black',anchor='nw'))
                            aligncount=aligncount+1
                            display_count=display_count+1
                            #print 'scanning, aligncount: ',aligncount,alignedseq[aligncount]
            #
            # Found the residue, now print
            #
            x=x_start+display_count*spacing

            colour_val=self.colour_function(calculation,aanum,restype)
            if colour_val=='nontit':
                colour='grey'
            elif colour_val:
                colour=self.get_colour_fromval(colour_val)
            else:
                colour='yellow'
                
            self.sequencesobjs.append(self.canv.create_rectangle(x,y,x+spacing,y+y_spacing,fill=colour,outline='white'))
            self.sequencesobjs.append(self.canv.create_text(x,y,text=aa1,fill='black',anchor='nw'))
            count=count+1
            display_count=display_count+1
            if alignedseq:
                aligncount=aligncount+1
        return

    #
    # ----
    #

    def get_colour_fromval(self,value):
        """Convert a value to a colour"""
        if abs(value)<=1.0:
            return 'green'
        elif value>2.0:
            return 'red'
        elif value>1.0:
            return 'pink'
        elif value<-2.0:
            return 'blue'
        elif value<-1.0:
            return 'cyan'

        return
