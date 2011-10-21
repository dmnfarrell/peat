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

from xml.sax import make_parser
from xml.sax.handler import ContentHandler 

class BSML_Handler(ContentHandler):
    
    def __init__(self):
        self.loading_seq=None
        return
    
    
    def startElement(self,name,attrs):
        if name=='Seq-data':
            self.loading_seq=1
            self.sequence=''
        return

    def endElement(self,name):
        if name=='Seq-data':
            self.loading_seq=None
        return
    
    def characters(self,ch):
        #
        # Record the chars if we should
        #
        import string
        if self.loading_seq:
            self.sequence=self.sequence+string.strip(ch)
        return
#
# -----
#

class sequence:
    #
    # This class is both used as stand-alone and as part of Protool
    # There might be some inconsistencies because of this
    #

    def __init__(self):
        return
    
    #
    # ---------
    #

    def read_DNA_sequence(self,filename):
        #
        # Read the sequence
        #
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        #for line in lines:
        #    print line,
        #
        # Figure out what kind of file it is - a bit limited at the moment
        #
        import string
        if lines[0][0]=='>' and string.upper(lines[0][1:3])!='DL':
            # FASTA
            #print 'Reading FASTA'
            return self.readfasta(filename)
        elif string.upper(lines[0][:3])=='>DL':
            # PIR
            #print 'READING pir'
            return self.readpir(filename)
        elif string.upper(lines[0][:5])=='LOCUS':
            # GenBank
            #print 'READING Genbank'
            return self.readGenBank(filename)
        elif lines[0][2:5]=='xml':
            for line in lines:
                if string.find(line,'DOCTYPE Bsml')!=-1:
                    return self.readBSML(filename)
            return None
        else:
            # FLAT file
            #print 'READING FLAT'
            return self.readflat(filename)

    #
    # ---------
    #        
        
    def readpir(self,filename):
        #
        # Read a pir file
        #
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        import string
        self.name=string.split(lines[0],';')[-1]
        self.sequence=string.strip(string.join(lines[2:],''))
        self.sequence=self.clean_seq(self.sequence)
        if self.sequence[-1]=='*':
            self.sequence=self.sequence[:-1]
        return self.sequence

    #
    # ------
    #

    def writepir(self,sequence,filename,name='Sequence'):
        """Save the sequence as a PIR file"""
        import os
        topdir=os.getcwd()
        fd=open(os.path.join(topdir,filename),'w')
        fd.write('>%s;\n\n' %name)
        buffer=''
        for char in sequence:
            buffer=buffer+char
            if len(buffer)==80:
                fd.write(buffer+'\n')
                buffer=''
        fd.write(buffer+'\n')
        fd.close()
        return
        
        
    #
    # ---------
    #

    def readfasta(self,filename):
        #
        # Read a fasta file
        #
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        import string
        self.name=string.strip(lines[0][1:])
        self.sequence=string.strip(string.join(lines[1:],''))
        self.sequence=self.clean_seq(self.sequence)
        if self.sequence[-1]=='*':
            self.sequence=self.sequence[:-1]
        return self.sequence
    #
    # ---------
    #

    def readflat(self,filename):
        #
        # Read a flat file
        #
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        import string
        self.sequence=string.strip(string.join(lines,''))
        self.sequence=self.clean_seq(self.sequence)
        return self.sequence    
    #
    # ---------
    #        
        
    def readGenBank(self,filename,lines=None):
        #
        # Read a GenBank file
        #
        if not lines:
            fd=open(filename)
            lines=fd.readlines()
            fd.close()
        count=0
        while lines[count][:6]!='ORIGIN':
            count=count+1
        #
        # DNA seq start
        #
        count=count+1
        import string
        self.sequence=string.join(lines[count:])
        self.sequence=self.clean_seq(self.sequence)
        self.sequence=string.replace(self.sequence,'//','')
        #
        # Get rid of all the numbers
        #
        for x in string.digits:
            self.sequence=string.replace(self.sequence,x,'')
        return self.sequence

    #
    # ---------
    #

    def readBSML(self,filename):
        #
        # Read a BSML file
        #
        parser = make_parser()   
        bsmlh = BSML_Handler()
        parser.setContentHandler(bsmlh)
        parser.parse(open(filename))
        return bsmlh.sequence
    #
    # ---------
    #

    def clean_seq(self,seq,remove_gaps=False):
        """
        # Do some standard cleaning
        """
        import string
        seq=string.strip(seq)
        seq=string.replace(seq,' ','')
        seq=string.replace(seq,'\n','')
        seq=string.replace(seq,'\r','')
        #
        # Odd chars?
        #
        seq=string.replace(seq,'\xbb','')
        seq=string.replace(seq,'\xbf','')
        seq=string.replace(seq,'\xef','')
        # 
        # Gaps?
        #
        if remove_gaps:
            seq=seq.replace('-','')
        seq=string.upper(string.strip(seq))
        return seq

    #
    # ---------
    #

    def length(self):
        return len(self.sequence)

    #
    # ---------
    #

    def count_hydrophobics(self):
        #
        # Returns the number of hydrophobics in the current molecule
        #
        # Hydrophobics are defined as: Ala, Phe, Ile, Leu, Met, Val, Tyr, Trp
        #
        hyd='AFILMVWY'
        count=0
        pir=self.PirSeq()
        for aa in pir:
            if aa in hyd:
                count=count+1
        return count
    #
    # ---------
    #

    def count_hydrophilic_and_charged(self):
        #
        # Returns the number of hydrophilic and charged amino acids in the current molecule
        #
        # Hydrophilic residues are: Ser, Thr, Asp, Asn, Gln, Glu, Lys, His, Arg
        #
        hyd='STDNQEKHR'
        count=0
        pir=self.PirSeq()
        for aa in pir:
            if aa in hyd:
                count=count+1
        return count
    #
    # ---------
    #

    def count_charged(self):
        #
        # Returns the number of charged residues in the current molecule
        #
        # Hydrophilic residues are: Asp, Glu, Lys, His, Arg
        #
        hyd='DEKHR'
        count=0
        pir=self.PirSeq()
        for aa in pir:
            if aa in hyd:
                count=count+1
        return count
    
    #
    # ---------
    #
    
    def PirSeq(self,seq=None):
        """
        # This function returns a string holding the sequence of the PDB file
        # (aa one-letter code)
        """
        import string
        pir=''
        if not seq:
            seq=self.sequence
        for residue,type in seq:
            if self.three_to_one.has_key(string.upper(type)):
                pir=pir+string.upper(self.three_to_one[string.upper(type)])
        return pir
    
    #
    # ---------
    #

    def Seq2Pir(self,seq=None,name='Dummy_protein',description=' '):
        """
        # Translates self.seq or seq into pir format.
        # The pir file is returned as a list of strings
        """
        import string
        pir=[]
        pir.append('>P1;'+name+'\n')
        pir.append(description+'\n')
        if not seq:
            seq=''
            for x in self.sequence:
                if self.three_to_one.has_key(x[1]):
                    seq=seq+self.three_to_one[x[1]]
        seq=string.upper(seq)
        seq=string.strip(seq)
        seq=string.replace(seq,' ','')
        seq=string.replace(seq,'\n','')
        x=0
        s=' '
        for char in seq:
            s=s+char
            x=x+1
            if x==60:
                pir.append(s+'\n')
                s=' '
                x=0
        pir.append(s+'*\n')
        return pir
    
    #
    # ---------
    #

    def Seq2FASTA(self,seq=None,ID='Dummy protein',description=' '):
        #
        # Return the sequence in fasta format
        #
        import string
        pir=self.Seq2Pir(seq,ID)
        fasta=['>'+ID+'|'+description]
        for line in pir[1:]:
            if string.strip(line)=='':
                pass
            else:
                line=string.replace(line,'*','')
                fasta.append(string.strip(line))
        return fasta

    #
    # ---------
    #

    def composition(self):
        #
        # Find the composition
        #
        comp={}
        for letter in self.sequence:
            if not comp.has_key(letter):
                comp[letter]=0
            comp[letter]=comp[letter]+1
        #
        letters=comp.keys()
        letters.sort()
        print 'Composition:'
        for letter in letters:
            print '%3s : %4d = %7.2f pc' %(letter,comp[letter],float(comp[letter])/float(self.length())*100.0)
        self.comp=comp
        return
    
    #
    # ---------
    #

    def subfragment(self,minlength=10,maxlength=None):
        #
        # Find if the sequence repeats itself
        #
        print 'Searching for exact repeats.'
        if not maxlength:
            maxlength=self.length()/2
        import string
        results={}
        for length in range(minlength,maxlength):
            for startpos in range(self.length()-length):
                fragment=self.sequence[startpos:startpos+length]
                if string.count(self.sequence,fragment)>1:
                    results[fragment]=string.count(self.sequence,fragment)
        print 'Substrings:'
        strings=results.keys()
        strings.sort()
        if len(strings)==0:
            print 'No exact repeats found.'
        for s in strings:
            print 'repeats: %3d : %s' %(results[s],s)
        return

    #
    # ------
    #
    
    def align_withseq(self,sequence):
        """Extract the sequence of the current molecule and align it with the sequence given.
        Return an array giving the relation between the sequence position and the PDB residue ID"""
        pdbpir=self.clean_seq(self.PirSeq(),remove_gaps=True)
        sequence=self.clean_seq(sequence,remove_gaps=True)
        import sequence_alignment
        NW=sequence_alignment.NW(pdbpir,sequence)
        al1,al2,map1,map2=NW.Align(verbose=True)
        map=[]
        residues=self.residues.keys()
        residues.sort()
        count=0
        #
        for count in range(len(sequence)):
            pdbrespos=map2[count]
            if pdbrespos!='-':
                map.append(residues[map2[count]])
            else:
                map.append(None)
        #
        # Check that we have 100% sequence identify for the aligned residues
        #
        print 'Aligned seqid',NW.aligned_seqid
        return map,NW.aligned_seqid

