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

"""needleman-Wunsch Algorithm for sequence alignment"""

import numpy as Numeric

class NW:
    """Needleman-Wunsch algorithm

    # The gap function, and the scoring functions need to be
    # enchanced to give a good alignment procedure.
    # Can be rewritten as Waterman-Smith relatively easily...
    """

    def __init__(self,seq1,seq2):
        self.s1=seq1
        self.s2=seq2
        return

    # 
    # -----
    #

    def NW_score(self,p1,p2,x_last=0,y_last=0):
        #
        # Simple scoring function
        #
        if None:
            score=0.0
            length_bonus=0.05
            while x_last!=0 and y_last!=0:
                if self.path[x_last,y_last]==3:
                    score=score+length_bonus
                    x_last=x_last-1
                    y_last=y_last-1
                else:
                    break
        #
        # Similarity in this position
        #
        score=0.0
        if self.s1[p1]==self.s2[p2]:
            score=score+1.0
        else:
            score=0.0
        return score

    #
    # ----
    #

    def NW_gap(self,x,y):
        #
        # Gap scoring function
        #
        if x==self.s1_len-1 or y==self.s2_len-1:
            return 0.0
        if self.path[x,y]==3:
            # Gap initiation
            return 1
        else:
            # Gap extension
            return 0.25

    #
    # -----
    #

    def Align(self, verbose=False):
        """
        # Needlenman-Wunsch alignment alg.
        """
        self.s1_len=len(self.s1)
        self.s2_len=len(self.s2)
        #
        # Set up the matrix
        #

        self.mat=Numeric.zeros((self.s1_len,self.s2_len))
        self.path=Numeric.zeros((self.s1_len,self.s2_len))
        #
        # Fill the edges
        #
        for x in range(self.s1_len):
            self.mat[x,0]=self.NW_score(x,0)
            self.path[x,0]=5
        for y in range(self.s2_len):
            self.mat[0,y]=self.NW_score(0,y)
            self.path[0,y]=6
        self.path[0,0]=0
        #
        # Fill the rest of the matrix
        #
        for x in range(self.s1_len)[1:]:
            for y in range(self.s2_len)[1:]:
                #
                # Score each position
                #
                up=self.mat[x,y-1]-self.NW_gap(x,y-1)+self.NW_score(x,y,x,y-1)
                left=self.mat[x-1,y]-self.NW_gap(x-1,y)+self.NW_score(x,y,x-1,y)
                diag=self.mat[x-1,y-1]+self.NW_score(x,y,x-1,y-1)
                if up>left:
                    if up>diag:
                        self.path[x,y]=2
                        self.mat[x,y]=up
                    else:
                        self.path[x,y]=3
                        self.mat[x,y]=diag
                else:
                    if left>diag:
                        self.path[x,y]=1
                        self.mat[x,y]=left
                    else:
                        self.path[x,y]=3
                        self.mat[x,y]=diag
        #
        # Done, now find the biggest edge value
        #
        x_pos=0
        y_pos=0
        max_score=-9999
        for x in range(self.s1_len):
            score=self.mat[x,self.s2_len-1]
            if score>max_score:
                x_pos=x
                y_pos=self.s2_len-1
                max_score=score
        for y in range(self.s2_len):
            score=self.mat[self.s1_len-1,y]
            if score>max_score:
                x_pos=self.s1_len-1
                y_pos=y
                max_score=score
        #
        # Found the end, build the overhang
        #
        al_s1=''
        ct_s1=[]
        al_s2=''
        ct_s2=[]
        if x_pos!=self.s1_len-1 and y_pos==self.s2_len-1:
            for x in range(self.s1_len-1,x_pos,-1):
                al_s1=self.s1[x]+al_s1
                ct_s1=['-']+ct_s1
                al_s2='-'+al_s2
        elif x_pos==self.s1_len-1 and y_pos!=self.s2_len-1:
            for y in range(self.s2_len-1,y_pos,-1):
                al_s2=self.s2[y]+al_s2
                ct_s2=['-']+ct_s2
                al_s1='-'+al_s1
        done=None
        last_x=0
        last_y=0
        while not done:
            #
            # Where did we come from?
            #
            step=self.path[x_pos,y_pos]
            #print step
            last_x=x_pos
            last_y=y_pos
            if step==1:
                x_pos=x_pos-1
                #print 'x-1'
            elif step==2:
                y_pos=y_pos-1
                #print 'y-1'
            elif step==3:
                x_pos=x_pos-1
                y_pos=y_pos-1
                #print 'diag'
            else:
                done=1
                break
            #
            # Put in the appropriate character
            #
            #print 'Adding'
            if last_x!=x_pos:
                al_s1=self.s1[last_x]+al_s1
            else:
                al_s1='-'+al_s1
            if last_y!=y_pos:
                al_s2=self.s2[last_y]+al_s2
            else:
                al_s2='-'+al_s2
            #
            # Update the ct arrays
            #
            if last_x!=x_pos and last_y!=y_pos:
                ct_s1=[last_y]+ct_s1
                ct_s2=[last_x]+ct_s2
            elif last_x!=x_pos and last_y==y_pos:
                ct_s1=['-']+ct_s1
            elif last_y!=y_pos and last_x==x_pos:
                ct_s2=['-']+ct_s2
            #print_pir_alignment(al_s1,al_s2) 
            #
            # If we reached the end for either of the sequences, then fill in the rest
            #
            if x_pos==0:
                #print 'filling in becuase x is 0'
                al_s1=self.s1[x_pos]+al_s1
                al_s2=self.s2[y_pos]+al_s2
                ct_s1=[y_pos]+ct_s1
                ct_s2=[x_pos]+ct_s2
                for y in range(y_pos-1,-1,-1):
                    al_s2=self.s2[y]+al_s2
                    ct_s2=['-']+ct_s2
                    al_s1='-'+al_s1
                break
            if y_pos==0:
                #print 'filling in becuase y is 0'
                al_s1=self.s1[x_pos]+al_s1
                al_s2=self.s2[y_pos]+al_s2
                ct_s1=[y_pos]+ct_s1
                ct_s2=[x_pos]+ct_s2
                for x in range(x_pos-1,-1,-1):
                    al_s1=self.s1[x]+al_s1
                    ct_s1=['-']+ct_s1
                    al_s2='-'+al_s2
                break

        if verbose == True:
            print '-------'
            print_pir_alignment(al_s1,al_s2)
        #
        # Calculate the sequence identity
        #
        identical=0
        id2=0
        totaligned=0 # Last two are counters measuring sequence identity of aligned positions
        for count in range(len(al_s1)):
            if al_s1[count]==al_s2[count]:
                identical=identical+1
            #
            if al_s1[count]!='-' and al_s2[count]!='-':
                if al_s1[count]==al_s2[count]:
                    id2=id2+1
                totaligned=totaligned+1
        #
        self.sequence_identity=float(identical)/len(al_s1)*100.0
        self.aligned_seqid=float(id2)/totaligned*100.0
        #
        # al_s2 and al_s2 are the aligned sequences
        # ct_s1 is a list of the aligned residues in s2. e.g. s2[ct_s1[x]] gives the s2 residue aligned to residue x in s1
        # ct_s2 is a list of the aligned residues in s1
        return al_s1,al_s2,ct_s1,ct_s2

#
# -----------
#

def pir2Protool(sequence):
    """ Reformats a regular one-letter sequence (e.g. 'ASDDE') to
        the protool type sequence ( [[':0001','ALA'], [':0002','SER'] ...)
    """
    import Protool, string
    X=Protool.structureIO()
    list=[]
    number=1
    for letter in sequence:
        list.append([':'+string.zfill(number,4),X.one_to_three[letter]])
        number=number+1
    return list

#
# ---------
#

def Protool2pir(sequence):
    """ Reformats a protool style sequence to a regular one-letter sequence"""
    import Protool
    X=Protool.structureIO()
    seq=''
    #
    # Get rid of final stop codon
    #
    if sequence[-1][1]=='***':
        sequence=sequence[:-1]
    ignored={}
    for aa in sequence:
        if aa[1]:
            #
            # Ignore water
            #
            if not X.three_to_one.has_key(aa[1]):
                ignored[aa[1]]=1
            else:
                seq=seq+X.three_to_one[aa[1]]
        else:
            seq=seq+'-'
    return seq,ignored

#
# -------
#

def print_pir_alignment(seq1,seq2):
    #
    # Print it
    #
    x=seq1   # new seq which will be aligned with the old one
    m=seq2   # old original seq
    #
    id=''
    for count in range(len(x)):
        if m[count]==x[count]:
            id=id+'*'
        else:
            id=id+' '
    line=80
    end=None
    count=0
    while not end:
        if (count+1)*line>=len(x):
            end=1
            print 'seq1',x[count*line:]
            print 'seq2',m[count*line:]
            print '    ',id[count*line:]
            print 
        else:
            print 'seq1',x[count*line:(count+1)*line]
            print 'seq2',m[count*line:(count+1)*line]
            print '    ',id[count*line:(count+1)*line]
            print
            count=count+1
    print
    return 
#
# -------
#

def read_clustal_aln_file(filename):
    """Read a clustalW/X aln file"""
    import os
    seqs={}
    if os.path.isfile(filename):
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        #
        # Skip first three lines
        #
        import string
        for line in lines[1:]:
            split=string.strip(line).split()
            if len(split)>1 and line[1]!=' ':
                if seqs.has_key(split[0]):
                    seqs[split[0]]=seqs[split[0]]+split[1]
                else:
                    seqs[split[0]]=split[1]
    return seqs
