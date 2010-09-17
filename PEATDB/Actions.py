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


"""
DB Actions class to handle peat actions from the within the GUI Application
Author: Damien Farrell 2010
"""

from Tkinter import *
import os
import tkSimpleDialog, tkFileDialog, tkMessageBox
from Base import zDatabase
from Record import PEATRecord
from Extfile import FileHandler

class DBActions(object):
    """This class handles misc application tasks like opening sub apps"""
    MUT = None
    yasara = None
    yasaraobjects = {}

    three_to_one={'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I',
                  'LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S',
                  'THR':'T','VAL':'V','TRP':'W','TYR':'Y','***':'*'}

    def __init__(self, DB=None):
       
        return


    @classmethod
    def AAList2String(self, aalist):
        """Convert the amino acid sequence as stored in aaseq to a string"""
        s = ''
        for a in aalist:
            s+=self.three_to_one[a[1]]
        return s
        
    @classmethod       
    def string2AAseq(self, seq, chain='A'):
        """Convert string to amino acid aaseq as stored in PEAT"""
        from PEAT_SA.Core import Utilities
        codemap  = Utilities.InvertedCodeDictionary()
        codes = [codemap[el] for el in list(seq)]
        indexes = [Utilities.CreateResidueCode(chain=chain, number=index) for index in range(1, len(seq) + 1)]        
        return zip(indexes, codes)
        
    @classmethod
    def mutationCodeFromSequences(self, initialSequence, targetSequence, offset=0, chain='A'):
        '''Returns a MutationSet representing the mutations necessary to change initialSequence to targetSequence            
        Parameters:
                initialSequence - "WildType"
                targetSequence - "Mutant"
                offset - The difference between the index of the first aa in the initial sequence and 1
                        e.g. if the first amino-acid in initialSequence string is actually the 5th
                        in the real sequence then the offset is 4.               
        Note: If this requires insertions or deletions this method returns None'''
        code=''
        muts=[]
        if len(initialSequence) != len(targetSequence):
            print 'Sequences are not the same length'
            return None
        if initialSequence.find('-') != -1:
            return None
        elif targetSequence.find('-') != -1:
            return None	
        for count in range(len(initialSequence)):
            initialAA=initialSequence[count]
            targetAA=targetSequence[count]
            if initialAA != targetAA:              
                muts.append(chain + str(count+1 + offset) + str(targetAA)) 
        code += '+'.join(muts)
        return code
    
    @classmethod
    def checkMutation(self, DB, name, ref=None):
        """Check mutations based on ref sequence and current mutant
           sequence, should be triggered whenever ref protein is altered so 
           that the mutation codes are updated.."""
        prot = DB.get(name)
        if prot.aaseq == None:
            return
        if ref == None:
            ref = self.DB.meta.refprotein
            
        refseq = self.AAList2String(DB.get(ref).aaseq)
        
        if prot.aaseq == None:
            return
        #get mutations from sequence
        seq = self.AAList2String(prot.aaseq)       
        if seq == refseq:
            return
        #assumes chain A..    
        mcode = self.mutationCodeFromSequences(refseq, seq, offset=0)
        prot.Mutations = mcode
        return
        
    def addProteinSeq(self, DB, name):
        """Add a protein sequence"""
        seq_win=Toplevel()
        seq_win.title('Enter amino acid sequence')
        lbl=Label(seq_win,text='Please enter AA sequence below (1-letter code) or click browse to select a file')
        lbl.grid(row=0,column=0,columnspan=8,sticky='w')
        seq_box=Text(seq_win,height=10,width=85)
        seq_box.grid(row=1,column=0,columnspan=15,sticky='news')
        self.seq_file=''
        self.seq_text=''
        self.seq_start=StringVar()
        self.seq_start.set("1")
        Label(seq_win,text='Sequence starts at number').grid(row=2,column=0,sticky='e',columnspan=2)
        Entry(seq_win,textvariable=self.seq_start,width=5).grid(row=2,column=2,sticky='w')

        def get_seq_file():
            import tkFileDialog
            savedir=self.preferences.get('datadir')
            self.seq_file=tkFileDialog.askopenfilename(defaultextension='.PIR',
                                                       initialdir=savedir,
                                                       filetypes=[("PIR file","*.PIR"),
                                                                  ("All files","*.*")],
                                                       parent=seq_win)

            import AA_sequence
            SEQ=AA_sequence.sequenceIO()
            SEQ.readpir(self.seq_file)
            seq_box.insert(END,SEQ.sequence)
            return
        Button(seq_win,text='Browse',command=get_seq_file).grid(row=2,column=4)
        def get_seq_text():
            self.seq_text=seq_box.get('1.0',END)
            seq_win.destroy()
        Button(seq_win,text='Create protein',command=get_seq_text).grid(row=2,column=6)
        def cancel_get_seq():
            self.seq_file=None
            seq_win.destroy()
            return
        Button(seq_win,text='Cancel',command=cancel_get_seq).grid(row=2,column=7)
        self.master.wait_window(seq_win)

        # Add the protein after checking

        if self.seq_file is None:
            #self.data['DBinstance'].delete_complete_record(protein_name)
            return

        import string
        if string.strip(self.seq_text)!='':
            import string
            sequence=string.strip(self.seq_text)
        else:
            import tkMessageBox
            tkMessageBox.showwarning('No sequence',
            'No sequence found\nNo protein added',
            parent=self.master)
            return

        # Change the sequence to three-letter code
        import DNAtool.mutation
        three_lt_seq=[]
        for aa in sequence:
            if DNAtool.mutation.one_to_three.has_key(aa.upper()):
                three_lt_seq.append(DNAtool.mutation.one_to_three[aa.upper()])
            else:
                import tkMessageBox
                tkMessageBox.showwarning('Invalid sequence',
                'Sequence contains invalid character: "%s"' %aa,
                parent=self.master)
                self.data['DBinstance'].delete_complete_record(protein_name)
                return

        # Add the protein
        ok=self.data['DBinstance'].add_protseq(protein_name,
                                               three_lt_seq,
                                               int(self.seq_start.get()))
        if not ok:
            raise 'Something went wrong when addding the protein sequence'
        return


    @classmethod
    def addPDBFile(self, DB=None, name=None, pdbfile=None, pdbdata=None):
        """Add a PDB file to the record given as argument"""
        import os        
        if pdbdata == None and pdbfile == None:
            savedir=os.getcwd()
            global PDB_code
            pdbfile=tkFileDialog.askopenfilename(defaultextension='.pdb',
                                         initialdir=savedir,
                                         filetypes=[("PDB file","*.pdb"),
                                                    ("PDB file","*.brk"),
                                                    ("All files","*.*")])
        #if not pdbfile or not pdbdata:
        #    return
        import Protool
        self.X=Protool.structureIO()
        # Extracting PDB_code from pdbfile
        if pdbdata != None: 
            self.X.readpdb(data=pdbdata)

        elif os.path.isfile(pdbfile):
            PDB_code=pdbfile.split('/').pop().split('.')[0]
            # Try to read it using Protool    
            try:
                self.X.readpdb(filename=pdbfile)
            except:
                import tkMessageBox
                tkMessageBox.showwarning('Error reading PDB file',
                                         'I could not read the PDB file. This probably means that the PDB file is corrupt in some way.')
                return

        # Extract the sequence
        import sequence_alignment
        pdb_1,ignored_res1=sequence_alignment.Protool2pir(self.X.sequence)
        print 'IGNORED',ignored_res1
        if ignored_res1!={}:
            igroups=ignored_res1.keys()
            igroups.sort()
            import tkMessageBox
            tkMessageBox.showwarning('Unknown entities in PDB file',
                                     'I ignored the following residue types/molecules in the PDB file:\n%s' %(str(igroups)))

        # Get the entry sequence
        accept_alignment_automatically=None
        record_AA = DB.get_AA_sequence(name)
        if record_AA:
            record_AA1,ignored_res=sequence_alignment.Protool2pir(record_AA)
        else:
            # If we do not have an amino acid sequence for the record, then
            # we simply use the one from the PDB file and accept the alignment
            # straight away
            accept_alignment_automatically=1
            import copy
            record_AA1=copy.deepcopy(pdb_1)

            # Also deposit the amino acid sequence in the protein record
            DB.data[name]['aaseq'] = copy.deepcopy(self.X.sequence)

        # Align the two sequences
        NW_align=sequence_alignment.NW(pdb_1,record_AA1)
        al_pdb,al_record,map_pdb,map_record=NW_align.Align()
        self.al_pdb = al_pdb
        self.al_record = al_record

        # Find regions of overlap

        ids=0
        for count in range(len(al_pdb)):
            res_pdb=al_pdb[count]
            res_rec=al_record[count]
            if res_pdb==res_rec:
                ids=ids+1
        print 'Sequence identity %5.3f' %(100.0*float(ids)/float(len(al_pdb)))
        self.AlignmentMap = {}
        self.AlignmentMap['OrigAa']=al_record
        self.AlignmentMap['AlignedAa']=al_pdb

        def store_PDB():
            DB.storePDB(name, self.X, self.AlignmentMap)
            AlignWindow.destroy()

        #Make alignment window
        AlignWindow=Toplevel()
        self.AlingWindow=AlignWindow
        AlignWindow.geometry('+100+200')
        AlignWindow.title('Please check alignment')
        AlignWindow.button = Button(AlignWindow,
                                    {"text": "Alignment OK", "fg": "black",
                                     "command":store_PDB})
        AlignWindow.button.grid(row=3,column=0)
        AlignWindow.button = Button(AlignWindow,
                                    {"text": "Alignment not OK", "fg": "black",
                                     "command": AlignWindow.destroy})
        AlignWindow.button.grid(row=3,column=1)
        AlignWindow.Slider=Scrollbar(AlignWindow,orient=HORIZONTAL)
        AlignWindow.Slider.grid(row=1,column=0,sticky='news',columnspan=2)

        listbox = Listbox(AlignWindow,{"height": 2,"width":80,"font":"courier 14"})
        listbox.insert('end',"PEAT_DB record: "+al_record)
        listbox.insert('end',"PDB file      : "+al_pdb)
        listbox.grid(row=0,column=0,columnspan=2)
        listbox.config(xscrollcommand=AlignWindow.Slider.set)
        AlignWindow.Slider.config(command=listbox.xview)

        if accept_alignment_automatically:
            store_PDB()
        return

    @classmethod    
    def fetchPDB(self, pdbid):
        import urllib
        url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid
       
        stream = urllib.urlopen(url).read()
        info = stream.info()
        status = info.status
        if status is not "":
            return None            
        elif not info.has_key('content-disposition'):
            return None
        elif len(pdb) < 10:
            return None
        else:
            return stream

    @classmethod
    def displayStructure(self, protein, field_name='Structure',
                         molapp=None, path=None, color='green', clear=True):
        """Display a structure"""

        pdblines,X = self.DB.getStructure(protein, field_name)
        if not pdblines:
            return
        pdbname=os.path.join(os.getcwd(),'._'+protein)
        # Write the pdbpipe
        fd=open(pdbname,'w')
        fd.writelines(pdblines)
        fd.close()
       
        if molapp == 'pymol':
            try:
                self.launchPyMol(pdbname)
            except:
                import subprocess
                print path,pdbname
                p=subprocess.Popen([path,pdbname])
        elif molapp == 'yasara':
            self.launchYasara(protein, pdbname, path, color, clear)

        elif molapp in ['vmd','rasmol']:
            import subprocess
            try:
                p=subprocess.Popen([molapp,pdbname])
            except:
                print 'command not found, trying path..'

                if path == '' or path == None:
                    print 'no path for application set...'
                else:
                    p=subprocess.Popen([path,pdbname])

        else:
            import tkMessageBox
            tkMessageBox.showwarning('Cannot display structure',
                                     'The structure cannot be displayed because no molecular '
                                     'graphics program has been set up.')
        return

    @classmethod
    def launchPyMol(self, pdbfile):
        """Open pdb in PyMOL"""
        import pymol
        # Call the function below before using any PyMOL modules.
        pymol.finish_launching()
        from pymol import cmd
        cmd.load(pdbfile)
        cmd.show_as('sticks')
        return

    @classmethod
    def launchYasara(self, protein, pdbfile, path, color='green', clear=False):
        """Launch yasara"""
        yasaradir = path
        import os
        if not os.path.isdir(yasaradir):
            yasaradir=os.path.split(yasaradir)[0]
        dirname=os.path.split(yasaradir)[1]
        if dirname.lower()=='yasara.app':
            yasaradir=os.path.join(yasaradir,'yasara')            
        else:
            pass         

        if self.yasara == None:
            self.importYasara(yasaradir)
        yasara = self.yasara

        # Get Yasara to load the PDB file
        if clear == True:
            yasara.Clear()
        yasara.run('OriAll 0,0,0')
        obj=yasara.LoadPDB(pdbfile)
        yasara.Style('stick')
        yasara.HideRes('Hoh')
        yasara.ColorObj(obj, color)
        yasara.ColorBG('black')
        self.yasaraobjects[pdbfile] = obj

        #If a reference protein is selected then try to colour the mutations
        DB = self.DB
        refprot = DB.meta.refprotein
        rec = DB.get(protein)
        if refprot is None:
            pass
        else:           
            parentrec = DB.get(refprot)
            is_parent, operations = rec.getAncestry(parentrec)
            if is_parent:
                for op in operations:
                    import pKa.pKD_tools
                    number=int(pKa.pKD_tools.get_resnum_from_mut(op))
                    yasara.ColorRes(number,'red')
           
        return 
    
    @classmethod
    def importYasara(self, yasaradir):
        import sys,os
        sys.path.append(os.path.join(yasaradir,'pym'))
        sys.path.append(os.path.join(yasaradir,'plg'))
        import yasaramodule as yasara
        self.yasara = yasara
        return
        
    @classmethod
    def initProtool(self, callback=None):
        """Import Protool"""
        # Init the modelling routines

        print 'Importing Protool'
        import Protool.mutate        
        self.MUT = Protool.mutate.Mutate()                    
        print 'Done importing Protool'
        return self.MUT
        
    @classmethod
    def writePDB(self, pdblines, filename):       
        fd=open(filename,'w')
        for line in pdblines:
            fd.write(line)
        fd.close()
        return
    
    @classmethod
    def makemutantSequence(self, sequence, operations):
        """Apply the specified mutations to a sequence and return the mutant seq
        Sequence must be in the [[A:0001:ALA],[A:0002:GLU]] format
        Operations is a list of the following types:
        Mutations: A:0001:ALA:ASP
        Deletions: delete:A:0002:GLU
        Insertions: insert:1:A:0003:THR:ALA, insert:2:A:0003:THR:TRP (insert THR,TRP after A:0003:THR)
        Operations are always performed in sequence numbering order  """
        
        if operations==[]:
            return sequence
        ops_sorted={}
        insertions=[]
        for operation in operations:
            s_op=operation.split(':')  
            # Normal mutation
            import pKa.pKD_tools as pKD_tools
            resid=pKD_tools.get_resid_from_mut(operation)
            if ops_sorted.has_key(resid):
                raise Exception('More than one operation on the same residue: %s' %resid)
            ops_sorted[resid]=['mutate',operation]
        
        # Perform the operations        
        new_seq=[]
        new_count=None
        new_chain=None
        for resid,restyp in sequence:
            # Make sure that the chain hasn't changed or if we are at the beginning then init
            if resid.split(':')[0]!=new_chain:
                #Initialise
                sp_resid=resid.split(':')
                new_chain=sp_resid[0]
                new_count=int(sp_resid[1])
                newresid='%s:%s' %(new_chain,string.zfill(new_count,4))
            # Does this residue have an operation?
            if ops_sorted.has_key(resid):
                op=ops_sorted[resid]
                if op[0]=='delete':
                    # Deletion
                    if op[1]==restyp:
                        pass # This deletes the residue
                    else:
                        raise Exception('Incorrect org residue in deletion: %s' %op)
                elif op[0]=='insert':
                    # Insertion
                    inserts=op[1].keys()
                    inserts.sort()
                    for i in inserts:
                        if i[0]==restyp:
                            new_seq.append([newresid,i[1]])
                            new_count=new_count+1
                            newresid='%s:%s' %(new_chain,string.zfill(new_count,4))
                elif op[0]=='mutate':
                    # Mutation
                    import pKa.pKD_tools as pKD_tools
                    orgres=pKD_tools.get_oldrestyp_from_mut(op[1])
                    if orgres==restyp:
                        new_seq.append([newresid,pKD_tools.get_newrestyp_from_mut(op[1])])
                        new_count=new_count+1
                        newresid='%s:%s' %(new_chain,string.zfill(new_count,4))
                    pass
                else:
                    raise Exception('Unknown mutations spec: %s' %op)
            else:
                new_seq.append([resid,restyp])
                new_count=new_count+1
                newresid='%s:%s' %(new_chain,string.zfill(new_count,4))
        return new_seq
    
    @classmethod
    def checkModels(self, DB=None, callback=None, selected=None, usemutationcodes=False):
        """Check that we have modelled a structure for everything we can"""
        if DB == None:
            return
        proteins = DB.getRecs()
        refprot = DB.meta.refprotein
        refseq = DB[refprot].aaseq
        refaa = self.AAList2String(refseq)
        refpdb = DB[refprot].Structure
        refpdbfile = os.path.join(os.getcwd(), 'ref.pdb')
        self.writePDB(refpdb, refpdbfile)
        failed = []

        # Check that Protool is loaded
        if not self.MUT:
            self.initProtool()
            
        #Create protool oinstance for ref pdb
        import Protool
        Xref = Protool.structureIO()
        Xref.parsepdb(refpdb)
            
        # Find all potential parents
        records_with_structure=[]
        for protein in proteins:
            rec = DB.get(protein)
            if rec.hasStructure() == 'available':
                records_with_structure.append(protein)
                
        # Loop over selected or all    
        if selected == None:
            selected = list(set(proteins) - set(records_with_structure))
        numrecords=len(selected)
        count=1
                    
        for protein in selected:
            rec = DB.get(protein)
            if rec.hasStructure() == 'available':
                continue                
            print 'Protein:', protein            
            
            #if no sequence try create one from mutation code
            if rec.aaseq == None and rec.Mutations != None:
                #print refaa
                print 'no sequence, using mutation code and ref protein seq'                 
                import PEAT_SA.Core as Core
                print 'Record has mutation code %s' %rec.Mutations
                mutationSet = Core.Data.MutationSet(rec.Mutations)
                mutseq = mutationSet.applyToSequence(refaa, id='A', pdb=Xref)                 
                rec.aaseq = self.string2AAseq(mutseq)
                     
            parent_with_structure = []
            for parent in records_with_structure:
                parentrec = DB.get(parent)                              
                is_parent, operations = rec.getAncestry(parentrec)
                         
                # We can only model on X-ray structures
                if parentrec.hasStructure() == 'available' and is_parent:                    
                    parent_with_structure.append([parent, len(operations)])

            # Record failure to model
            if parent_with_structure == []:
                continue

            # Find the best parent
            def compare_func(x,y):
                if x[1]>y[1]:
                    return 1
                elif x[1]==y[1]:
                    return 0
                if x[1]<y[1]:
                    return -1
                
            parent_with_structure.sort(cmp=compare_func)
            parent = parent_with_structure[0][0]               
            operations = rec.getAncestry(parentrec)[1]
            print 'Using %s as template with %d operations.' %(parent, len(operations))

            # Start the modelling
            pdblines = parentrec.Structure
            # Load the pdb file
            import Protool
            X=Protool.structureIO()
            X.parsepdb(pdblines)
            self.MUT.new_PDB(X)
            self.MUT.max_tolerated_bump=0.5
            atom_changes=[]
            skip_protein=None
            self.MUT.new_mutation()
            
            for operation in operations:                
                # Is this a deletion?                    
                if operation.find('delete')!=-1:
                    print 'This is a deletion - Jens should write code for modelling this'
                    print 'Deletion ignored for now'
                    continue
                elif operation.find('insert')!=-1:
                    print 'This is an insertion - Jens should write code for modelling insertions'
                    print 'Insertion ignored for now'
                    continue
                
                # This is a normal mutation                    
                # Get the residue number, old residue and new residue
                
                import pKa.pKD_tools as pKD_tools
                new_res = pKD_tools.get_newrestyp_from_mut(operation)
                old_res = pKD_tools.get_oldrestyp_from_mut(operation)
                resid = pKD_tools.get_resid_from_mut(operation)
                
                #import string
                if not X.residues.has_key(resid):
                    print 'No structural info for mutation %8s. Not modelling this mutation\n' %operation
                    continue
                
                # Actually make the mutation                
                bump_score=self.MUT.mutate(resid,new_res,orgtype=old_res)
                print 'Mutation: %s, bump_score: %s' %(resid+new_res,str(bump_score))
                if bump_score is None:
                    skip_protein=True
                    break
                else:
                    atom_changes=atom_changes+self.MUT.mutate_operations
                self.MUT.mutate_operations=[]

            # Update progress            
            completion = float(count)/float(numrecords)*100.0
            if callback != None:
                callback(completion)
            else:
                print '%4d of %4d, completion; %5.2f%%' %(count,float(numrecords),completion)
            count=count+1
            
            # Did it work?                
            if skip_protein:         
                print
                print 'Modelling failed for %s' %protein
                failed.append(protein)
                rec.Structure = 'Bumps'
                rec.structuretype = 'failed model'
                continue
            
            # We have all sets of changes in atom_changes           
            rec.Structure = {'Rotamer_operations': atom_changes}
            rec.Structure['parent'] = parent
            rec.structuretype = 'peat model'
            
        print 'Done'
        if len(failed)>0:
            print 'Failed to model the following proteins:'
            for f in failed: print f
        return

    @classmethod
    def modelFromMutationCode(self):
        """Model directly from mutation code"""
        import PEAT_SA.Core as Core
        print 'Record has mutation code %s' %rec.Mutations
        mutationSet = Core.Data.MutationSet(rec.Mutations)
        mutationCodes = mutationSet.mutationCodes(Xref, reduced=False)
        #print refpdbfile, mutationCodes
        result, score = Protool.mutate.Model_Mutations(refpdbfile, [],
                        mutationCodes,return_score=True)
        if result == False:
            print 'failed to model'
        else: 
            mutant = result.PI
            mutant.writepdb('mutant.pdb')                
            #rec.Structure = mutant.write_pdb('dummy',nowrite=1)
            rec.structuretype = 'protool model'
        return
        
    @classmethod
    def getRecordsSelector(self, DB, dwin):
        """Selection boxes for record/col selection"""
        recs = DB.getRecs()
        fields = DB.getFields()

        yscrollbar=Scrollbar(dwin,orient='vertical',width=12)
        yscrollbar.grid(row=1,column=2,sticky='news',padx=2)
        recsbox=Listbox(dwin,bg='white',
                         exportselection=0,
                         height=18,
                         width=20,
                         yscrollcommand=yscrollbar.set,
                         selectmode=EXTENDED)
        yscrollbar.config(command=recsbox.yview)
        Label(dwin,text='Records:').grid(row=0,column=0,sticky='news')
        recsbox.grid(row=1,column=0,columnspan=2,sticky='news',padx=1,pady=3)
        recsbox.config(state=NORMAL)

        y1scrollbar=Scrollbar(dwin,orient='vertical',width=12)
        y1scrollbar.grid(row=1,column=5,sticky='news',padx=2)
        colsbox=Listbox(dwin,bg='white',
                         exportselection=0,
                         height=20,
                         width=20,
                         yscrollcommand=y1scrollbar.set,
                         selectmode=EXTENDED)

        y1scrollbar.config(command=colsbox.yview)
        Label(dwin,text='Fields:').grid(row=0,column=3,sticky='news')
        colsbox.grid(row=1,column=3,columnspan=2,sticky='news',padx=1,pady=3)
        colsbox.config(state=NORMAL)
        dwin.rowconfigure(1, weight=1)
        dwin.columnconfigure(0, weight=1)
        for r in recs:
            recsbox.insert('end', r)
        for f in fields:
            colsbox.insert('end', f)
        return recsbox, colsbox            
