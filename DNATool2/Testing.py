#!/usr/bin/env python
#
# DNATool - A program for DNA sequence manipulation
# Copyright (C) 2012- Damien Farrell & Jens Erik Nielsen
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
# Email: farrell.damien_at_gmail.com 


"""Test functions for DNATool"""

import os, sys, types
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord

def writeFastaFile(sequences):
    recs=[]
    for s in sequences:
        recs.append(SeqRecord(Seq(s)))
    SeqIO.write(recs, "test.faa", "fasta")
    return

def clustalAlignment(filename):
    from Bio.Align.Applications import ClustalwCommandline
    cline = ClustalwCommandline("clustalw", infile=filename)
    print 'performing alignment..'    
    stdout, stderr = cline()
    align = AlignIO.read("test.aln", "clustal")
    '''print align
    from Bio import Phylo
    tree = Phylo.read("test.dnd", "newick")
    Phylo.draw_ascii(tree)'''
    return align

clustalAlignment('test.faa')


