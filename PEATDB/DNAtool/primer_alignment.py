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

"""Functions for dealing with the primer database in DNAtool"""

from Tkinter import *
import string
import mutation

class align_primer:

    """Functions for aligning primers to template DNA sequences"""

    def get_real_primer_position(self,this_primer):
        """Get the real base position for a primer"""
        rev_compl=None
        position=this_primer['startpos']
        primer_sequence=this_primer['sequence']
        if not primer_sequence or primer_sequence=='':
            return
        if position<0:
            rev_compl=1
            position=len(self.parent.data['DNAseq'])-(abs(position)+len(primer_sequence))
            #
            # Reverse
            #
            primer_sequence=''
            for letter in this_primer['sequence']:
                primer_sequence=letter+primer_sequence
            #position=position+1
        return position

    #
    # This draws a primer OR any comparison sequence above the main sequence
    # Use tags to distinguish primers and other kinds of sequences

    def draw_primer(self,this_primer,colour='blue',thetag='primer',seqname=None,lvl=0):
        """This draws a primer OR any comparison sequence above the main sequence"""

        # Draw the primer
        font = self.parent.getCurrentFont()
        objs={}
        mismatches=0
        rev_compl=None
        position=this_primer['startpos']
        primer_sequence=this_primer['sequence']
        if not primer_sequence or primer_sequence=='':
            return
        if position<0:
            rev_compl=1
            position=len(self.parent.data['DNAseq'])-(abs(position)+len(primer_sequence))

            # Reverse
            primer_sequence=''
            for letter in this_primer['sequence']:
                primer_sequence=letter+primer_sequence
        position=position+1
        DNA_stretch=''
        match=''

        # Write the 5' or 3'
        x,y=self.parent.get_base_pos_on_screen(position-2)

        # First find the y value, corrected for scale and no. of sequences etc.
        # This needs to be only done once, so use the dummy ytmp value when doing
        # get_base_pos_on_screen after this

        y=(y-10)/self.parent.y_scale
        #move primer further up, if the comparison DNA sequence is displayed
        if self.parent.show_comp_sequence.get()==1 and thetag=='primer':
            y=(y-10)/self.parent.y_scale
            #raise up primer correct no. of 'levels' os it is just above highest one
            if lvl>0:
                lvl=lvl-1
                y=y-(lvl*15)/self.parent.y_scale
        #raise mutliple levels when many comparison sequences displayed
        if thetag != 'primer':
            y=y-(lvl*15)/self.parent.y_scale

        if this_primer['startpos']<0:
            start_symbol="3'"
            end_symbol="5'"
        else:
            start_symbol="5'"
            end_symbol="3'"

        objs[self.parent.seqframe.create_text(x,y,text=start_symbol,font=font,
                                            fill='blue',anchor='w',tag=thetag)]=1

        # Now draw the primer/seq
        for letter in primer_sequence:
            x,ytmp=self.parent.get_base_pos_on_screen(position)
            # If we are still within the parent sequence then detect differences

            if position-1<len(self.parent.data['DNAseq']):
                org_base=self.parent.data['DNAseq'][position-1]
            else:
                org_base=letter

            identical=None
            if rev_compl:

                if mutation.match(org_base,letter):
                    identical=1
            else:
                if string.upper(org_base)==string.upper(letter):
                    identical=1
            #
            # Draw primer/seq. If the seq has a name then add the tag for that.
            #
            if identical:
                if seqname==None:
                    obj=objs[self.parent.seqframe.create_text(x,y,text=letter,font=font,
                                                        fill=colour,anchor='w',
                                                        tags=thetag)]=1
                else:
                    name=seqname.rstrip('seq.clipped')
                    obj=objs[self.parent.seqframe.create_text(x,y,text=letter,font=font,
                                                        fill=colour,anchor='w',
                                                        tags=(thetag,name))]=1
                match=match+'|'
            else:
                if seqname==None:
                    obj=objs[self.parent.seqframe.create_text(x,y,text=letter,font=font,
                                                        fill='red',anchor='w',
                                                        tags=thetag)]=1
                else:
                    name=seqname.rstrip('seq.clipped')
                    obj=objs[self.parent.seqframe.create_text(x,y,text=letter,font=font,
                                                        fill='red',anchor='w',
                                                        tags=(thetag,name))]=1

                match=match+' '
                mismatches=mismatches+1
            DNA_stretch=DNA_stretch+org_base
            position=position+1

        # Write the end symbol
        x,ytmp=self.parent.get_base_pos_on_screen(position)
        objs[self.parent.seqframe.create_text(x,y,text=end_symbol,font=font,
                                                fill='blue',anchor='w',tags=thetag)]=1
        return mismatches,match,objs

    def find_primer_binding_sites(self,primer_seq,DNA_seq):

        '''Find all possible binding sites for the primer in the DNA_seq
         - also probe the reverse strand'''

        binding_sites=self.find_primer_binding_sites_1(primer_seq,DNA_seq)
        # Reverse complement strand
        rev_compl = mutation.get_reverse_complementary(DNA_seq)
        rev_binding_sites = self.find_primer_binding_sites_1(primer_seq,rev_compl,1)
        return binding_sites+rev_binding_sites


    def find_primer_binding_sites_1(self,primer_seq,DNA_seq,rev=None):
        """Find all possible binding sites for the primer in the DNA_seq.
        This function should be changed so we use the scan_primer class in mutation.py"""
        #

        scores=[]
        primer_len=len(primer_seq)
        for pos in range(0,len(DNA_seq)):
            DNA_stretch=DNA_seq[pos:min(pos+primer_len,len(DNA_seq))]
            score=0
            for comp_pos in xrange(min(len(primer_seq),len(DNA_stretch))):
                if DNA_stretch[comp_pos]==primer_seq[comp_pos]:
                    score=score+1
            scores.append([score,pos])
        #
        # Find the ten best binding sites
        #
        scores.sort()
        vals=scores[-10:]
        vals.reverse()
        best_ten=[]
        for pos in vals:
            if not rev:
                best_ten.append([pos[1],pos[0]])
            else:
                best_ten.append([-pos[1],pos[0]])
        return best_ten[:3]
    #
    # -----------
    #

class dummy_align_primer(align_primer):

    def __init__(self,parent):
        self.parent=parent
        return

