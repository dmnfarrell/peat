#!/usr/bin/env python
#
# pKa - various programs and scripts for pKa value analysis, calculation and redesign
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

def read_aln(msf_file):
    import msf_parser
    seqs=msf_parser.read_msf(msf_file)
    #
    print 'I have read: %s ' %msf_file
    seq_keys=seqs.keys()
    seq_keys.sort()
    print 'The file contains %d sequences.' %(len(seq_keys))
    #
    # Strip all gaps
    #
    #    import string
    #    nogap={}
    #    for seq in seqs.keys():
    #        nogap[seq]=string.replace(seqs[seq],'.','')
    return seq_keys,seqs

def intro_text():
    #
    # Print intro text
    #
    print
    print '     pKa_alignment'
    print '-----------------------'
    print
    print 'Create pka setup for a large number of PDB files from a multiple sequence alignment'
    print
    return

#
# ----------------------------------
#

def usage_text():
    print
    print 'Usage: pKa_alignment.py <msf_file> [<pdb1] <pdb2> .. <pdbn>]'
    print 'All sequences IDs in the msf file must correspond to a pdb file'
    print 'in the pwd'
    print 'If <msf_file> is -, then the sequences of all pdb files given are'
    print 'extracted and dumped into a big pir file'
    print
    print
    return

#
# ----------------------------------
#

def clean_files(pdbfiles):
    #
    # Clean all PDB files
    #
    dict={}
    print
    print 'Cleaning all PDB files (deleting water molecules and reconstructing missing side chains)'
    print
    import WI_tools, os
    for pdb in pdbfiles:
        print 'Cleaning: %s' %pdb
        newfile=pdb+'.clean.pdb'
        dict[pdb]=newfile
        if not os.path.isfile(pdb+'.clean.pdb'):
            logfile,files=WI_tools.corall_delwat(os.path.join(os.getcwd(),pdb))
            if len(files.keys())!=1:
                raise "Could not clean ",pdb
            #
            # Write the new pdbfile
            #
            
            fd=open(newfile,'w')
            for line in files[files.keys()[0]]:
                fd.write(line)
            fd.close()
    print
    print 'All files cleaned'
    print '----------------------------------------'
    return dict
#
# ----------------------------------
#

def read_sequences(pdbfiles,newfiles):
    #
    # Load all the PDB files and Make sure that the sequences are ok
    #
    align={}
    allseqs={}
    print 'Reading pdb files and extracting sequences'
    import Protool
    for pdb in pdbfiles:
        pdbfile=newfiles[pdb]
        X=Protool.structureIO_fast()
        print 'Reading: %s' %pdbfile
        X.readpdb(pdbfile)
        X.RemoveALT()
        s_keys=X.residues.keys()
        s_keys.sort()
        #
        # Construct a special list for the sequence
        #
        newlist=[]
        for s in s_keys:
            if X.three_to_one.has_key(X.resname(s)):
                newlist.append([s,X.three_to_one[X.resname(s)]])
        align[pdb]=newlist[:]
        allseqs[pdb]=X.Seq2Pir(None,pdb)
    return align,allseqs

#
# ------------------------------------------
#
        
def main():
    intro_text()
    #
    # Get the filename
    #
    import sys
    try:
        msf_file=sys.argv[1]
    except:
        usage_text()
        raise Exception()
    #
    # if the filename starts with - then we only want ro extract all the sequences
    #
    makeseqs=None
    if msf_file[0]=='-':
        #
        # We only want to extract sequences from the pdb files
        #
        # The rest of the argv contains the pdb files
        #
        pdbfiles=sys.argv[2:]
        newfiles=clean_files(pdbfiles)
        align,allseqs=read_sequences(pdbfiles,newfiles)
        #
        # Write pir file
        #
        seqs=allseqs.keys()
        seqs.sort()
        import os
        pirfile=os.path.join(os.getcwd(),'all_seqs.pir')
        fd=open(pirfile,'w')
        for seq in seqs:
            for line in allseqs[seq]:
                fd.write(line)
            fd.write('\n')
        fd.close()
        #
        # Done
        # Now the user has to align that file
        #
        return
    #
    # ----------------------------------------------
    #
    else:
        #
        # This the normal case. Setup the pKa calculations
        #
        #
        # Read the alignment
        #
        seq_keys,seqs=read_aln(msf_file)
        newfiles=clean_files(seq_keys)
        align,allseqs=read_sequences(seq_keys,newfiles)
        #
        # Check that the sequences correspond to the alignment sequences
        # AND build an alignment of the pdb sequences
        #
        pdb_align={}
        for seq in align.keys():
            pdbcount=0
            pdb_align[seq]=[]
            for letter in seqs[seq]:
                if letter=='.':
                    pdb_align[seq].append('.')
                else:
                    if letter!=align[seq][pdbcount][1]:
                        print
                        print 'Different letter at position %d ' %pdbcount
                        print 'in %s ' %seq
                        print 'alignment says: %1s while pdbfile says: %1s' %(letter,align[seq][pdbcount][1])
                        print
                        print seqs[seq]
                        print align[seq]
                        raise "sequences do not match"
                    else:
                        print align[seq][pdbcount]
                        pdb_align[seq].append(align[seq][pdbcount][0])
                        pdbcount=pdbcount+1
        print 'All sequences ok'
        print
        #
        # Prompt the user to get the file to work with
        #
        print 'Select the one that you want to identify active site groups for'
        for seq in seq_keys:
            print seq
        print
        done=None
        while not done:
            define_seq=raw_input('Enter name: ')
            if seqs.has_key(define_seq):
                done=1
        #
        # ok use selected a sequence
        #
        # open the corresponding pdbfile and print the sequence in there
        #
        for residue in align[define_seq]:
            print residue[0],residue[1]
        print
        input=raw_input("Enter residues (separate by ,): ")
        import string
        residues=string.split(input,',')
        #
        # Check that the residues exist
        #
        align_numbers=[]
        for residue in residues:
            found=None
            number=0
            for res in pdb_align[define_seq]:
                #print residue,res
                if residue==res:
                    found=1
                    align_numbers.append(number)
                    break
                number=number+1
            if not found:
                raise "Residue not found: ",residue
        print 'All residues found in %s' %define_seq
        print align_numbers
        for number in align_numbers:
            print number,pdb_align[define_seq][number]
        #
        # ok, prepare all the directories
        #
        for seq in align.keys():
            import os
            dir=seq+'_pka'
            if not os.path.isdir(dir):
                os.mkdir(dir)
            #
            # Copy the pdbfile
            #
            pdbfile=os.path.join(os.getcwd(),dir,seq+'.clean.pdb')
            if os.path.isfile(pdbfile):
                os.unlink(pdbfile)
            os.link(os.path.join(os.getcwd(),seq+'.clean.pdb'),pdbfile)
            #
            # Find the number of the titratable group
            #
            import pKa
            Y=pKa.pKatools()
            groups=Y.get_titratable_residues2(seq+'.clean.pdb')
            groupnums=[]
            #
            # Get the residues from the alignment
            #
            residues=[]
            for number in align_numbers:
                residues.append(pdb_align[seq][number])
            print 'I found these residues %s for %s' %(str(residues),seq)
            for residue in residues:
                found=None
                for group in groups:
                    print group[0],residue
                    print
                    print 'You have to fix this for the new output from get_titratable_residues'
                    print
                    raise Exception()
                    # !!!!!!
                    if group[0]==residue:
                        groupnums.append(group[2])
                        found=1
                        break
                if not found:
                    import Protool
                    T=Protool.structureIO_fast()
                    T.readpdb(seq+'.clean.pdb')
                    print 'I could not find a titratable group to match this residue: %s %s in %s' %(residue,T.resname(residue),seq)
                    a=raw_input('Should I just ignore this group? (Y/N) ')
                    if a=='y' or a=='Y':
                        pass
                    else:
                        raise 'Could not identify group'
            #
            # Write the Invocation file
            #
            grps=''
            for grp in groupnums:
                grps=grps+str(grp)+','
            grps=grps[:-1]
            invocation=os.path.join(os.getcwd(),dir,'Invocation')
            fd=open(invocation,'w')
            fd.write('/net/home/jnielsen/lib/python/pKarun.py %s -dbcrit 1000 -subset -group_cutoff 2.0 -groups %s\n' %(seq+'.clean.pdb',grps))
            fd.close()
            #
            # Copy the parameter files
            #
            os.system('cp /net/home/jnielsen/pkaparms/TOPOLOGY.H '+dir)
            os.system('cp /net/home/jnielsen/pkaparms/DELRAD.DAT '+dir)
            os.system('cp /net/home/jnielsen/pkaparms/DELCRG.DAT '+dir)
        return
    

if __name__=="__main__":
    main()
