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
from seq_utils import *
import PEATDB.DNA_sequence as DNA_sequence
import re, sys, os, string

# conversion from three to one-letter code

three_to_one={'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I',
              'LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S',
              'THR':'T','VAL':'V','TRP':'W','TYR':'Y','***':'*'}
#
# Conversion from one-letter to three-letter code
#
one_to_three={}
for tlt in three_to_one.keys():
    one_to_three[three_to_one[tlt]]=tlt
#
# Genetic code to amino acid dictionary
#
genetic_code={'TTT':'PHE', 'TTC':'PHE','TTA':'LEU','TTG':'LEU','TCT':'SER','TCC':'SER','TCA':'SER','TCG':'SER',
              'TAT':'TYR', 'TAC':'TYR','TAA':'***','TAG':'***','TGT':'CYS','TGC':'CYS','TGA':'***','TGG':'TRP',
              'CTT':'LEU','CTC':'LEU','CTA':'LEU','CTG':'LEU','CCT':'PRO','CCC':'PRO','CCA':'PRO','CCG':'PRO',
              'CAT':'HIS','CAC':'HIS','CAA':'GLN','CAG':'GLN','CGT':'ARG','CGC':'ARG','CGA':'ARG','CGG':'ARG',
              'ATT':'ILE','ATC':'ILE','ATA':'ILE','ATG':'MET','ACT':'THR','ACC':'THR','ACA':'THR','ACG':'THR',
              'AAT':'ASN','AAC':'ASN','AAA':'LYS','AAG':'LYS','AGT':'SER','AGC':'SER','AGA':'ARG','AGG':'ARG',
              'GTT':'VAL','GTC':'VAL','GTA':'VAL','GTG':'VAL','GCT':'ALA','GCC':'ALA','GCA':'ALA','GCG':'ALA',
              'GAT':'ASP','GAC':'ASP','GAA':'GLU','GAG':'GLU','GGT':'GLY','GGC':'GLY','GGA':'GLY','GGG':'GLY'}
#
# -----------------------------
#

def get_codons(AA):
    """Get all the codons that code for amino acid AA"""
    
    posskey = []
    for codon in genetic_code.keys():   
        if genetic_code[codon].upper()==AA.upper():
            posskey.append(codon)
    return posskey

#
# -----------------------------
#

#
# Inverse genetic code
#
inverse_genetic_code={}
for key in three_to_one.keys():
    inverse_genetic_code[key]=get_codons(key)
#
# Define DNA complementarity
#
complementary={'A':'T','T':'A','C':'G','G':'C'}


try:
    from Numeric import *
except:
    from numpy import *


class restriction_digest:
    """Class for performing restriction digests
       Uses precompiled regular expressions to get maximum speed"""

    def __init__(self):
        """Initialise"""
        #
        # Find the file with the restriction enzymes
        #
        found=None
        for scriptdir in sys.path:
            restriction_enzyme_file=os.path.join(scriptdir,'restriction_enzymes.DAT')
            if os.path.isfile(restriction_enzyme_file):
                found=1
                break
            restriction_enzyme_file=os.path.join(scriptdir,'DNAtool','restriction_enzymes.DAT')
            if os.path.isfile(restriction_enzyme_file):
                found=1
                break
        #
        # Did we find the file?
        #
        if not found:
            print
            print 'File holding description of restriction enzymes could not be found'
            print
            return 
        #
        # Open the file
        #
        fd=open(restriction_enzyme_file)
        lines=fd.readlines()
        fd.close()
        #
        # Parse the lines
        #
        self.restriction_lines=[]
        #
        # Is this a REBASE file?
        #
        found=None
        count=0
        while not found:
            header=lines[count].strip().split()
            if len(header)>0:
                found=1
            count=count+1
        if header[0]!='REBASE':
            import tkMessageBox
            tkMessageBox.showwarning('Not a REBASE file',
                                     'The restriction enzyme definition file is not in REBASE format')
            return
        #
        # OK, REBASE file - this is the only format we deal with right now
        #
        self.parse_REBASE(lines)
        #
        # Vars
        #
        self.exclude_promiscuous_enzymes=1
        #
        # Pre-compile all the regular expressions
        #
        self.compile_regexs()
        #
        # Success
        #
        return

    #
    # ----------
    #

    def parse_REBASE(self,lines):
        """Parse a REBASE restriction enzyme file"""
        count=1
        while 1:
            while (lines[count][:3]!='<1>' or lines[count][3]=='<') and count+1<len(lines):
                count=count+1
            if not count+1<len(lines):
                return
            #
            # At start of record
            #
            name=lines[count][3:].strip()
            rec_seq=lines[count+2][3:].strip()
            self.restriction_lines.append('%s %s' %(name,rec_seq))
            count=count+1
        return
            

    #
    # -----------
    #

    def compile_regexs(self):
        #
        # Compile all the regular expressions and store them in a ditionary (self.enzymes_regexs)
        #
        self.enzymes_regexs={}
        import re, string
        #
        # Loop over the enzyme definitions
        #
        for line in self.restriction_lines:
            #
            # Get the enzyme name
            #
            enz_name=line.split()[0]
            #
            # Compile the regex
            #
            rec_seq=line.split()[1]
            #
            # Deal with the (x/y) statements
            # These mean that we add max(x,y) Ns to the sequence and insert cuts at those positions
            #
            while rec_seq.find('(')!=-1:
                p_start=rec_seq.find('(')
                p_end=rec_seq.find(')')
                txt=rec_seq[p_start:p_end+1]
                nums_t=txt[1:-1].split('/')
                nums=[]
                for char in nums_t:
                    nums.append(int(char))
                max_N=max(nums)
                second_site=min(nums)
                new_txt='^'+max_N*'N'
                if p_start!=0:
                    rev=''
                    for l in new_txt:
                        rev=l+rev
                    new_txt=rev
                    new_txt=new_txt[:second_site]+'^'+new_txt[second_site:]
                else:
                    new_txt=new_txt[:max_N-second_site+1]+'^'+new_txt[max_N-second_site+1:]
                rec_seq=rec_seq.replace(txt,new_txt)
            #
            # Find all positions where the enzyme cuts
            #
            poscuts=[]
            count=0
            for l in rec_seq:
                if l=='^':
                    poscuts.append(count-len(poscuts))
                count=count+1
            #
            # If no position is specified, then we guess at the middle
            #
            if poscuts==[]:
                poscuts=[int(len(rec_seq)/2)]
            #
            # Produce the regular expression
            #
            rec_seq =  re.sub ( r"\^", r"", rec_seq )
            rec_seq=self.expand_recognition_sequence(rec_seq)
            #
            # We only inlude enzymes with N's in their recognition sequence if the flags are right
            # 
            if rec_seq.find('N')>0:                
                if self.exclude_promiscuous_enzymes==1:
                    continue
            comp_seq=re.compile(rec_seq)
            #
            # Add the enzyme
            #
            self.enzymes_regexs[enz_name]={'re':comp_seq,'realseq':rec_seq,'poscut':poscuts,
                'ncuts':len(poscuts)}
        #
        # Done
        #
        return 

        #
        # ----------
        #
        
    def expand_recognition_sequence(self,seq):
        #
        # Expand all the aliases
        #
        old_seq=seq
        seq = string.replace(seq,'R','[AG]')
        seq = string.replace(seq,'Y','[TC]')
        seq = string.replace(seq,'W','[AT]')
        seq = string.replace(seq,'S','[GC]')
        seq = string.replace(seq,'M','[AC]')
        seq = string.replace(seq,'K','[GT]')
        #
        seq = string.replace(seq,'B','[GCT]')
        seq = string.replace(seq,'H','[ATC]')
        seq = string.replace(seq,'V','[GAC]')
        seq = string.replace(seq,'D','[GAT]')
        #
        seq=seq.replace('N','[GATC]')
        #
        # If the flag is set and we have a change in sequence
        # then we shouldn't use this enzyme
        #
        if self.exclude_promiscuous_enzymes==1:
            if old_seq!=seq:
                return 'Do not include'
        return seq
        
    #
    # ---------------------
    #

    def get_restriction_sites(self,DNA_sequence,enzyme_list):
        #
        # Get all the restriction sites
        #
        if enzyme_list!=None:
            #
            # If we get a certain list of enzymes then use that one
            #
            sites=self._get_restriction_sites(DNA_sequence,enzyme_list)
        else:
            #
            # Just use all enzymes
            #
            sites=self._get_restriction_sites(DNA_sequence,self.list_of_enzymes)
        return sites

    #
    # -----------
    #

    def _get_restriction_sites(self,DNA_seq,enzyme_list):
        #
        # Find the restriction sites by matching regular expressions
        #
        sites={}
        for enzyme_name in enzyme_list:
            if self.enzymes_regexs.has_key(enzyme_name):
                reg_ex=self.enzymes_regexs[enzyme_name]['re']
                #
                # Find the sites
                #
                iterator=reg_ex.finditer(DNA_seq)
                #
                # Loop over each site (each site is a re.match object in iterator)
                #
                for match in iterator:
                    #
                    # do we have an entry for this enzyme already?
                    #
                    if not sites.has_key(enzyme_name):
                        sites[enzyme_name]=[]
                    #
                    # Find the exact position of the cut and append it to the list
                    #
                    for cut in self.enzymes_regexs[enzyme_name]['poscut']:
                        sites[enzyme_name].append([cut+match.start()])
            #
            # done
            #
        return sites


# =============================================================================


def match(base1,base2):
    #
    # Check the bases pair
    #
    if base1==complementary[base2]:
        return 1
        
    return None

#
# --------
#

def get_reverse_complementary(DNA_seq):
    #
    # Produce the reverse complementary DNA sequence
    #
    inverted_DNA_compl=''
    for base in DNA_seq:
        inverted_DNA_compl=complementary[base]+inverted_DNA_compl 
    return inverted_DNA_compl

#
# ---------
#
        
class scan_primer:

    """Scan the DNA sequence for binding sites for a primer"""

    def __init__(self,DNA_seq):
        #
        # Store the DNA sequence and the inverted complementary DNA sequence
        #
        import string
        self.DNA_seq=string.upper(DNA_seq)
        #
        # Make the inverted complementary sequence
        #
        self.inverted_DNA_compl=get_reverse_complementary(DNA_seq)
        #
        # done
        #
        return

    #
    # --------
    #

    def match(self,base1,base2):
        #
        # Check the bases pair
        #
        if match(base1,base2):
            return 1
        return -1

    #
    # ----------
    #

    def scan(self,primer):
        #
        # Scan over the DNA sequence to see if we have any other good binding sites
        #
        bad_pos = {}
        bad_pos_DNA_seq = {}
        for pos in range (0,(len(self.inverted_DNA_compl)-len(primer))):
            selfscorematch = 0
            DNA_seq_selfscorematch = 0
            for j in range (0, len(primer)):
                selfscorematch=selfscorematch+self.match(primer[j],self.inverted_DNA_compl[pos+j])
                DNA_seq_selfscorematch=DNA_seq_selfscorematch+self.match(primer[j],self.DNA_seq [pos+j])
            #
            # Done looping for this starting position
            #
            sum_inverted_DNA_compl=selfscorematch
            sum_DNA_seq=DNA_seq_selfscorematch
            #
            # Store the resulst in a dictionary
            #
            if sum_inverted_DNA_compl >= 8:
                bad_pos [pos] = {primer:sum_inverted_DNA_compl}
            if sum_DNA_seq >= 8:
                bad_pos_DNA_seq [pos] = {primer:sum_DNA_seq}
        return bad_pos, bad_pos_DNA_seq 

#
# ==================================================================
#

def initialise():
    #
    # Get the restriction enzyme definitions and the test sequence
    #
    s=getRNA_seq()
    return s,enzymes


def getRNA_seq():
    #
    # Read the RandomRNA.seq file
    #
    fd = open ('testDNA', 'r')
    DNAseqs = fd.readlines ()#Reads lines and returns the lines in a list
    fd.close ()
    import string
    s = string.join (DNAseqs,'')
    s=string.upper(s)
    return s

#
# -------------------
# 


def check_DNA(DNA_sequence):
    """Check that we have a DNA sequence without junk"""
    #
    import string
    # Remove all spaces
    DNA_sequence=string.replace(DNA_sequence,' ','')
    # Upper case
    DNA_sequence=string.upper(DNA_sequence)
    # Check that we only have DNA bases in the seq
    ok=1
    garbage={}
    DNA_bases=['A','G','C','T']
    for letter in DNA_sequence:
        if not letter in DNA_bases:
            ok=None
            garbage[letter]=1
    if ok:
        return ok, DNA_sequence
    return ok,garbage.keys()

#
# -------------------
# 

def translate(DNA_sequence):
    """Translate the DNA sequence in three forward frames"""
    #
    # Check the sequence first
    #
    ok,DNA_sequence=check_DNA(DNA_sequence)
    if not ok:
        print
        print 'DNA sequence not ok'
        print DNA_sequence
        print
        return None
    #
    # Translate
    #
    import string
    translation_3=[]
    translation_1=[]
    for start in range(3):
        position=start
        thisAA=''
        thisframe3=[]
        thisframe1=''
        while position<len(DNA_sequence):
            thisAA=thisAA+DNA_sequence[position]
            position=position+1
            if len(thisAA)==3:
                thisframe1=thisframe1+three_to_one[string.upper(genetic_code[thisAA])]
                thisframe3.append(string.upper(genetic_code[thisAA]))
                thisAA=''
        translation_3.append(thisframe3)
        translation_1.append(thisframe1)
    return translation_3,translation_1



#
# ----------------
#
    

def PCR_primer(DNA_seq,Tm_desired=65):
    SP = 1
    FP = 100
    s = DNA_seq
    #
    #megaprimer design
    #
    megaprimer_dict = {}
    for changedbase in [SP,FP]:
        #
        # Search in 3' durection (forward)
        #
        for i in range (0, len(s)):
            if changedbase-i>=0:
                primer_backward = s[changedbase-i:changedbase]
                AUcount, GCcount, Tm_backward, base_count_backward = AT_GC_Tm_base_count (primer_backward)                
            if changedbase+i<=len(s):
                primer_forward = s[changedbase:changedbase+i]
                AUcount, GCcount, Tm_forward, base_count_forward = AT_GC_Tm_base_count (primer_forward)
                if (Tm_backward+Tm_forward >= Tm_desired) and (base_count_forward + base_count_backward) % 3 == 0:
                    Tm = Tm_backward+Tm_forward
                    primer = primer_backward + primer_forward
                    megaprimer_dict [changedbase] = primer
                    break        
    #    
    #inverts the primer starting in FP
    #
    inverted = ''
    for count in range(len(megaprimer_dict[FP])):
        pos=-count-1
        inverted=inverted+megaprimer_dict [FP][pos]
    megaprimer_dict [FP] = inverted
    #
    # Tjuhej
    #
    return megaprimer_dict

#
# --------
#

def get_primer_Tm(DNAseq,primer,primer_start,method='Stratagene'):
    """Get the melting temperature of a primer. The entire primer must
    be contained in the DNA sequence"""
    #
    # If the starting position <0, then it means that the primer binds to the complementary strand
    #
    if primer_start<0 and primer_start!=None:
        DNAseq=get_reverse_complementary(DNAseq)
        primer_start=abs(primer_start)
    #
    # Precalculate the characteristics
    #
    if len(primer)==0:
        return 0.0,0.0,'None'
    if DNAseq:
        #
        # We can do this if we have an alignment
        #
        mismatch=0
        GC_count=0
        DNA_pos=primer_start
        for base in primer:
            temp_base=DNAseq[DNA_pos]
            if base!=temp_base:
                mismatch=mismatch+1
            if temp_base=='G' or temp_base=='C':
                GC_count=GC_count+1
            DNA_pos=DNA_pos+1
    else:
        #
        # No alignment
        #
        #GC_count=get_base_count(primer, 'GC')
        
        for base in primer:
            if base=='G' or base=='C':
                GC_count=GC_count+1
        
    perc_GC=float(GC_count)/float(len(primer))*100.0
    #
    # Different Tm calculation formulas
    #
    if method=='Stratagene' and DNAseq!=None and primer_start!=None:
        #
        # Calculate Tm using: Tm = 81.5 + 0.41(%GC) - 675/N - %mismatch
        #
        perc_mismatch=float(mismatch)/float(len(primer))*100.0
        return 81.5+0.41*perc_GC-675.0/float(len(primer))-perc_mismatch,mismatch,'Stratagene'
    elif method=='Simple' or DNAseq==None:
        #
        # Calculate Tm using: Tm =   81.5 + 0.41(%GC) - 675/N
        #
        return 81.5+0.41*perc_GC-675.0/float(len(primer)),0.0,'Simple'

    #elif method =='Alt' or DNASeq==None:
        #
        # test method: Tm = 81.5 + 16.6log[Na+] + 0.41(%G+C) - 600/N
        #
    elif method =='Basic' or DNASeq==None:
        #
        # Use the simplest, Tm = 2AT + 4GC
        # 
        AT=get_base_count(primer, 'AT')
        GC=get_base_count(primer, 'GC')
        return 2*AT + 4*GC,0.0,'Basic'        
        
    return 0.0,0.0,'None'


#
# ---------
#

def hairpin_propensity (primer):
    #
    # score for hairpin formation ###
    # according to Vallone and Butler: BioTechniques 37:226-231 (August 2004)
    #
    # Rewritten completely by Jens to give correct results 24/8-05
    #
    # Hairpin needs at least four bases for looping.
    # Count the number of matches on either side of the loop
    #
    scoredict = {}
    #
    # Test with loop lengths from 4 to 5
    #
    for loop_length in range(4,6):
        #
        # Try different positions for the loop
        #
        for pos in range (2, (len(primer)-(loop_length-1))):
            #
            # Get the DNA sequence before and after the loop
            #
            b4loop = primer[:pos]
            b4loop_r=''
            for char in b4loop:
                b4loop_r=char+b4loop_r
            #
            afterloop = primer[pos+loop_length:]
            length=min(len(b4loop),len(afterloop))
            score = 0    
            for i in range(length):
                base_before=b4loop_r[i]
                base_after=afterloop[i]
                if match(base_before,base_after):
                    score=score+1
                else:
                    score=score-1
            scoredict[str(pos)+':'+str(loop_length)] =  score
    #
    # Find the worst score
    #
    scorelist=scoredict.values()
    #
    # Make sure that we can deal with primers that are very short
    #
    maxscorelist=0
    if len(scorelist)>0:
        maxscorelist = max (scorelist)
    return scoredict,maxscorelist

#
# ----------
#

    
def self_complementarity_propensity(primer):
    #
    # score for self-complementary sequence ###
    # according to Vallone and Butler: BioTechniques 37:226-231 (August 2004)
    #
    # This function needs to be tested to see if it is accurate
    #
    
    selfscoredict = {}
    inverted_primer = ''
    partial_primer = ''
    matrix=zeros((len(primer)*len(primer)))
    score_dict = {}
    score=[]
    matrix.shape=(len(primer),len(primer))
    for count in range(len(primer)):
        pos=-count-1
        inverted_primer=inverted_primer+primer[pos]
    ##construction of the matching matrix##
    for j in range (0,len(primer)):
        for i in range (0,len(inverted_primer)): 
           if (primer [j] == 'A' and inverted_primer [i] == 'T') or (primer [j] == 'C' and inverted_primer [i] == 'G') or (primer [j] == 'G' and inverted_primer [i] == 'C') or (primer[j] == 'T' and inverted_primer[i] == 'A'):
                 matrix[j,i]= 1
           if (primer [j] == 'A' and inverted_primer [i] != 'T') or (primer [j] == 'C' and inverted_primer [i] != 'G') or (primer [j] == 'G' and inverted_primer [i] != 'C') or (primer[j] == 'T' and inverted_primer[i] != 'A'):
                 matrix[j,i]= -1
                    
                   
    for j in range(-len(primer),len(primer)):
        score_dict[j]=(trace(matrix,j))
        score.append(trace(matrix,j))
    #
    #

    score.append(0) # Prevent errors when we have an empty sequence []
    sum = max(score)
    return score_dict, sum

#
# ----------
#

def exhaustive_research(DNA_sequence,AA_number,NewAA,Tm_desired,Tm_method,enzyme_list=None,
                        parent=None,primer_max_len=30,find_restriction_site=None):
    """Generate a primer for introducing the required mutation by generating all possible DNA_sequences
    that will produce the amino acid sequence, and then scoring all of them for Tm, number of mismatches and
    restriction sites introduced"""
    #
    # Make the mutation in the DNA
    #
    codons=get_codons(NewAA)
    base_position=(AA_number-1)*3
    mut_DNA=DNA_sequence[:base_position]+'XXX'+DNA_sequence[base_position+3:]
    #
    # Get the Amino acid sequence from 5 aa (15 bases) up and downstream if possible
    #
    DNA_start=base_position-int(0.5*primer_max_len)
    DNA_end=base_position+int(0.5*primer_max_len)
    if DNA_start<0:
        DNA_end=DNA_end+abs(DNA_start)
        DNA_start=0
    if DNA_end>len(DNA_sequence):
        DNA_end=len(DNA_sequence)
    #
    # Find the base count of the mutation
    #
    primer_template=mut_DNA[DNA_start:DNA_end]
    mut_pos=primer_template.find('XXX')
    primer_template=primer_template.replace('XXX',codons[0])
    AA_sequence3,AA_sequence_1=translate(primer_template)
    AA_sequence3=AA_sequence3[0] # We are only interested in the first reading frame
    #
    # Keep a record of the template DNA
    #
    DNA_template=DNA_sequence[DNA_start:DNA_end]
    #
    # Generate all base sequences that will produce this amino acid sequuence
    #
    # We start by introducing as few mutations as possible, and end of with ?
    #
    res_comb = []
    max_permutations=20000
    #for codons_in_play in range(1,11):
    #    pass
    #print 'Mutated AA',AA_sequence3[AA_number]
    #print AA_sequence
    #stop
    for item in AA_sequence3:
        res_comb.append(inverse_genetic_code[item])   
    #
    # Call add_next
    #
    primers_comb=add_next(res_comb)
    #
    # Now we search for all the restriction sites on all the primers. Compare them with the ones present in the
    # DNAsequence and creates two different Dictionaries for them.
    #
    characteristics={}
    #
    # Instantiate the scan_primer class
    #
    SCAN=scan_primer(DNA_sequence)
    #
    # Loop over all primers
    #
    print 'Found %d primers in total' %(len(primers_comb))
    count=0
    parent.prog_can.create_text(5,5,text="Diversity generation",anchor='nw')
    box=None
    for comb in primers_comb:
        count=count+1
        #
        # The progress bar
        #
        if parent.cancelled:
            return None,None,None
        parent.update_primer_progress(float(count)/float(len(primers_comb)),1)
        #
        # Find sub-fragments of this primer that have the correct Tm
        #
        primer_length=len(comb)
        min_primer_length=15
        for primer_start in range(0,12,3):
            for primer_end in range(1,13,3):
                primer_seq=comb[primer_start:-primer_end]
                Tm,mismatches,Tm_method_used=get_primer_Tm(DNA_template,primer_seq,primer_start,Tm_method)
                if abs(Tm-Tm_desired)<4.0:
                    #
                    # Scan for other binding sites
                    #
                    characteristics[primer_seq]={'Tm':Tm,'mismatches':mismatches}
                if Tm<Tm_desired-4.0:
                    break
            #
            if primer_end<3 and Tm<Tm_desired-4.0:
                break
    print 'I found %7d combinations with a Tm less than 4K from the desired Tm' %(len(characteristics.keys()))
    print '--------------------------------------'
    #
    # Perform the resitriction site analysis
    #
    # Instantiate the restriction site class and get the restriction sites for the wt sequence
    #
    RS=restriction_digest() # This one calls the __init__ function
    rest_sites_DNATOT=RS.get_restriction_sites (DNA_sequence, enzyme_list)
    #
    # Finds all the restriction enzymes for each primer
    #
    rest_enzymes_primer = {}
    for primer in characteristics.keys():
        rest_enzymes_primer[primer]=RS.get_restriction_sites(primer,enzyme_list)
    #
    # compare the cut positions in the primer and in the DNAsequence
    #
    enzymes_that_already_cut={}
    new_enzymes={}
    #
    # Print the text
    #
    parent.prog_can.create_text(5,45,text="Primer restriction analysis",anchor='nw')
    count=0
    #
    #
    #
    for primer in characteristics.keys():
        #
        # The progress bar
        #
        if parent.cancelled:
            return None,None,None
        count=count+1
        parent.update_primer_progress(float(count)/float(len(characteristics.keys())),2)
        #
        # Find the enzymes that cut
        #
        for key_primer_enzymes in rest_enzymes_primer[primer].keys():
            if rest_sites_DNATOT.has_key(key_primer_enzymes):
                if not enzymes_that_already_cut.has_key (primer):
                    enzymes_that_already_cut [primer] = []
                enzymes_that_already_cut [primer].append (key_primer_enzymes)
            else:
                if not new_enzymes.has_key (primer):
                    new_enzymes [primer] = []
                new_enzymes[primer].append (key_primer_enzymes)
    #
    # Primer characterisation
    #
    primers_results_dict = {}
    counter=0
    #
    # Calculate self complementarity and hairpin-forming propensity
    #
    parent.prog_can.create_text(5,85,text="Primer characterisation",anchor='nw')
    count=0
    for primer in characteristics.keys():
        #
        # Progress bar
        #
        if parent.cancelled:
            return None,None,None
        count=count+1
        parent.update_primer_progress(float(count)/float(len(characteristics.keys())),3)
        #
        # Get the characteristics
        #
        hairpin_dict, maxscorelist = hairpin_propensity (primer)
        complementarity_dict,maxselfscorelist = self_complementarity_propensity (primer)
        if maxselfscorelist < 8:
            maxselfscorelist = 'OK' 
        elif maxselfscorelist >= 8:
            maxselfscorelist = 'BAD'
        if maxscorelist < 8:
            maxscorelist = 'OK'
        elif maxscorelist >= 8:
            maxscorelist = 'BAD'
        #
        # Add the characteristics to the dictionary
        #
        characteristics[primer]['hairpin_prop']=maxscorelist
        characteristics[primer]['self-compl']=maxselfscorelist
    return characteristics, new_enzymes, enzymes_that_already_cut


#
# ---------
#
           
    
def add_next (res_comb):
    """ Generates all the possible permutations by calling itself"""
    l = []
    for item in res_comb [0]:
        if len (res_comb) > 1:
            new = add_next (res_comb [1:])
            for item2 in new:
                l.append (item + item2)
        else:
            l.append (item)
    return l


#
# --------
#

def usage():
    print
    print 'Usage: mutation.py -seq <DNA sequence file> -mutation <mutation>'
    print 'mutation is of form D87N, where D is the org AA, 87 is the AA number, and N is the new AA'
    print 'The first letter of the DNA sequence file must be the first base of the starting codon'
    print '-fasta: DNA sequence is in fasta format'
    print '-pir: DNA sequence is in pir format [default]'
    print '-Tm <value>: Desired Tm for the primer'
    print '-no_restriction_site: Do not design a silent mutation that will introduce or eliminate a restriction site'
    print
    return

#
# -------
#

def get_defaults():
    #
    # Types for parameter passing:
    # text: text. res_num_list: List of residu numbers + associated number for each
    # number: number
    # res_list: List of residue numbers
    # T/F: True/False: If you give the argument then set to true, otherwise false
    # 
    defaults={'seq':['HEWL.fasta','text'],'mutation':['L10C','text'],'pir':[1,'T/F'],
              'fasta':[None,'T/F'],'Tm':[65.0,'number'],'no_restriction_site':[None,'T/F']}
    return defaults.copy()

#
# -------
#

def main():
    #
    # Parse the arguments
    #
    import sys,string
    defaults=get_defaults()
    #
    args=string.join(sys.argv[1:])
    args=string.split(args,'-')
    for arg in args:
        split=string.split(string.strip(arg))
        if split==[]:
            continue
        parm_name=split[0]
        if not defaults.has_key(parm_name):
            raise 'Unknown parameter: ',parm_name
        #
        # Deal with T/F
        #
        if len(split)==1:
            if defaults[parm_name][1]=='T/F':
                if defaults[parm_name][0]:
                    defaults[parm_name][0]=None
                else:
                    defaults[parm_name][0]=1
        #
        # Deal with all the other cases
        #
        elif len(split)==2:
            if defaults[parm_name][1]=='number':
                defaults[parm_name][0]=string.atof(split[1])
            else:
                defaults[parm_name][0]=split[1]
        else:
            raise 'Incorrect usage'
    #
    # Reformat
    #
    params={}
    for key in defaults.keys():
        params[key]=defaults[key][0]
    #
    # Load the file and design the primer
    #
    if not params['seq']:
        usage()
    seq_file=params['seq']
    new_AA=params['mutation'][-1]
    AA_number=int(params['mutation'][1:-1])
    
    
    S=DNA_sequence.sequence()
    DNA_seq=S.readpir(seq_file)
    Tm_desired=params['Tm']
    find_restriction_site=1
    if params['no_restriction_site']:
        find_restriction_site=None
    #
    # Do it!
    #
    import Protool
    X=Protool.structureIO()
    new_AA=X.threetoone[new_AA]
    #
    # Call the function for designing primers
    #
    new_enzymes, primers_results_dict, enzymes_that_already_cut, primer_starting_position,comb_on_Tm = exhaustive_research(DNA_seq,AA_number,new_AA,Tm_desired,find_restriction_site,enzyme_list=None)
    megaprimer_dict =  megaprimer (DNA_seq,Tm_desired=65)
    return 

#
# --------
#
   

if __name__=='__main__':
    #import profile
    #profile.run('main()','DNA')
    #import pstats
    #p = pstats.Stats('DNA')
    #p.strip_dirs().sort_stats(-1).print_stats()
    #p.sort_stats('cumulative').print_stats(20)
    #try:
    main()
    #except:
    #    usage()


    
