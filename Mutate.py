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
Class to handle creation of on the fly point mutants in PEAT
Author: Jens Nielsen, rewritten by Damien Farrell 2010
"""

from Tkinter import *
import os
import tkSimpleDialog, tkFileDialog, tkMessageBox
from Actions import DBActions

class PointMutate(object):
    """This class provides dialogs and functions for creating mutants"""
    def __init__(self, DB=None):
        self.DB = DB
        return

    def addMutant(self, DB, name, parentframe=None):
        """Open the 'Add mutant protein' dialog"""

        def close_mut_win(event=None):
            self.enter_mut=None
            self.mut_win.destroy()
            return

        self.DB = DB
        if parentframe!=None:
            self.mut_win = Frame(master=parentframe,relief=RAISED)
            self.mut_win.pack(fill=BOTH)
        else:
            self.mut_win=Toplevel()
            self.mut_win.title('Add mutant protein')
            self.mut_win.geometry("+200+200")
   
        row=0
        lbl=Label(self.mut_win,text='Parent protein:')
        lbl.grid(row=row,column=0,pady=5)
        self.parent = StringVar()

        self.parent_button=Menubutton(self.mut_win,textvariable=self.parent,relief=RAISED)
        self.parent_button.grid(row=row,column=1,sticky='news',pady=5)
        self.parent_menu=Menu(self.parent_button,tearoff=0)
        self.parent_button['menu']=self.parent_menu
      
        proteins = list(DB.getRecs())
        proteins.sort()
        self.parent.set(proteins[0])
        for protein in proteins:
            self.parent_menu.add_radiobutton(label=DB[protein]['name'],
                                             variable=self.parent,
                                             value=DB[protein]['name'],
                                             indicatoron=0,
                                             command=self.update_mut_selection)
        
        # Find the max and minimum residue number if we have a parent        
        maxres=1
        minres=1
        if self.parent.get()!='' and self.parent.get()!='Select parent':
            protein=self.parent.get()
            self.old_parent = protein
          
            # Can we get info from DNAtool?            
            if DB[protein].has_key('ORF_selected'):
                minres = DB[protein]['ORF_selected']['aastart_number']
                maxres = minres + DB[protein]['ORF_selected']['length']
            elif DB[protein].has_key('aaseq'):
                pass
            else:
                aaseq=[]
            
            # Dig out the Chain IDs
            
            #if self.data['DBinstance'].DB[protein].has_key('Structure'):
            #    pdblines,X=self.get_structure(protein,'Structure')
            #    import Protool
            #    X=Protool.structureIO()
            #    X.parsepdb(pdblines)
          
            sequence=DB.get_protein_sequence(protein)
            if  sequence:
                aaseq=sequence
                CIDs={}
                for resnum,aa in aaseq:
                    x=resnum.split(':')
                    if not CIDs.has_key(x[0]):
                        CIDs[x[0]]=1
                CIDs=CIDs.keys()
                CIDs.sort()
            else:
                CIDs=['']
        else:
            CIDs=['']
        
        # Get the widgets for specifying muts from DNAtool        
        row=row+2
        import DNAtool
        self.DT=DNAtool.MainWindow(data_passed=None,
                                   openwin=None)
        row=self.DT.specify_mut_res(window=self.mut_win,
                                    row=row,
                                    command=self.update_mut_selection,
                                    minres=minres,
                                    maxres=maxres,
                                    ChainIDs=CIDs)
        
        # Field for the name        
        row=row+1
        l=Label(self.mut_win,text="Enter name")
        l.grid(row=row,column=0, sticky='news', padx=2, pady=2)
        
        # Suggest a name for the mutant        
        name=self.suggest_name(self.parent.get(),
                                          self.DT.resnum.get(),
                                          self.DT.new_res.get(),
                                          self.DT.resname.get(),
                                          self.DT.chainID.get())
        self.new_name=StringVar()
        self.new_name.set(name)
        e=Entry(self.mut_win, textvariable=self.new_name, state=NORMAL, bg='white')
        e.grid(row=row,column=1,sticky='news', padx=2, pady=2)
        
        # Submit button        
        row=row+1
        b=Button(self.mut_win,text='Cancel',command=close_mut_win)
        b.grid(row=row,column=1,sticky='news', padx=2, pady=2)
        b2=Button(self.mut_win,text='Add mutant',command=self.enter_mut_data)
        b2.grid(row=row,column=0,sticky='news', padx=2, pady=2)
        b2=Button(self.mut_win,text='Add mutants from file',command=self.enter_file_muts)
        b2.grid(row=row+1,column=0,sticky='news', padx=2, pady=2)
        
        # Checkbox for showing the structure        
        self.show_struct=StringVar()
        self.show_struct.set(0)
        self.zoom_on_res=IntVar()
        self.zoom_on_res.set(1)
        Label(self.mut_win,text='Zoom on Residue').grid(row=row+2,column=0, padx=2, pady=2)
        Checkbutton(self.mut_win,onvalue=1,offvalue=0,command=self.update_mut_selection,
                      variable=self.zoom_on_res).grid(row=row+2,column=1)
        '''if DBActions.yasara:
            Label(self.mut_win,text='Show structure').grid(row=row+2,column=2, padx=2, pady=2)
            Checkbutton(self.mut_win,onvalue=1,offvalue=0,command=self.update_mut_selection,
                        variable=self.show_struct).grid(row=row+2,column=3)
            self.show_struct.set(1)
        else:
            Label(self.mut_win,text='Start Yasara').grid(row=row+3,column=0, padx=2, pady=2)
            Checkbutton(self.mut_win,onvalue=1,offvalue=0,command=self.open_yasara,
                          variable=self.show_struct).grid(row=row+3,column=1)'''
        import sys
        sys.stdout.flush()
        
        # Bindings        
        self.mut_win.bind('<Return>',self.update_mut_selection)
        #self.mut_win.bind('<Key>',self.update_mut_selection)
        self.DT.mutscale.configure(command=self.update_mut_selection)
        self.DT.mutscale.configure(from_=1)
        self.DT.mutscale.configure(to=1)
        
        # Wait for something        
        #self.master.wait_window(self.mut_win)
        
        # self.enter_mut is not created on error        
        if not getattr(self,'enter_mut',None):
            return
        if self.enter_mut:
   
            # Adding the ':' here is a dirty hack and I should do it properly by using these residue numbers
            # in DNAtool (jesus)            
            
            import string
            resnumber=':'+string.zfill(self.DT.resnum.get(),4)
            newmutation='%s:%s:%s:%s' %(self.DT.chainID.get(),
                                      resnumber.replace(':',''),
                                      self.DT.resname.get(),
                                      self.DT.new_res.get())
            new_name=self.new_name.get()
            self.data['DBinstance'].file_mutation_data(parent_name=self.parent.get(),mutation=newmutation,protein_name=new_name)
            #self.update_display()
            self.update_table()
        return


    def enter_file_muts(self,event=None):
        """Enter mutations from a file"""
        import tkMessageBox
        tkMessageBox.showwarning('Enter mutations from file',
                                'Note that the format of the mutations must be like'
                                'the following: '
                                'Single mutation: MUT1 A:0020:GLY:GLU'
                                'Multiple mutations: MUT2_DBL A:0020:GLY:GLU+A:0021:ASP:ALA'
                                'One mutant protein per line.',
                                parent=self.mut_win)

        import tkFileDialog
        savedir=self.preferences.get('datadir')
        filename=tkFileDialog.askopenfilename(defaultextension='.*',
                                    initialdir=savedir,
                                    filetypes=[("All files","*.*")],
                                    parent=self.mut_win)
        if not filename:
            return
        
        # Open the file        
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        for line in lines:
            import string
            line=string.strip(line)
            name=line.split()[0]
            mutations=line.split()[1]
            if self.data['DBinstance'].DB.has_key(name):
                tkMessageBox.showwarning('Enter mutations from file',
                                        'A protein with the name: "%s" already exists.\nSkipping this mutant.',
                                        parent=self.mut_win)
            else:
                self.data['DBinstance'].file_mutation_data(parent_name=self.parent.get(),mutation=mutations,protein_name=name)
        self.mut_win.destroy()
        return

    def enter_mut_data(self,event=None):
        """Process the data from the 'Add mutant' dialog and enter it into the Database"""

        orgres=None
        aa_number=self.DT.resnum.get() # Get the residue number
        seq=self.DB[self.parent.get()].aaseq
        if seq:
            for res in seq:
                num=int(res[0].split(':')[1])
                if num==aa_number:
                    orgres=res[1]
        else:
            import tkMessageBox
            tkMessageBox.showwarning('No sequence',
                                     'This parent protein does not have a sequence.\nI cannot add this mutant',
                                     parent=self.mut_win)
            return
        if not orgres:
            seq.sort()
            import tkMessageBox
            tkMessageBox.showwarning('Invalid sequence number',
                                     'Please give a sequence number from %s to %s' %(seq[0][0].split(':')[1],
                                                                                     seq[-1][0].split(':')[1]))
            return
        
        # Check the name        
        protein_name=self.new_name.get()
        if self.data['DBinstance'].is_known_protein(protein_name):
            import tkMessageBox
            tkMessageBox.showwarning('Name already in use',
                                     '%s is already in use. Please select a unique name' %protein_name)
            return 
        
        # Enter the data        
        self.enter_mut=1
        self.mut_win.destroy()
        return
    

    def suggest_name(self,parent,aa_number,newresidue,org_res,ChainID):
        
        """Suggest a name for new mutation"""
        
        if not parent:
            return '?+?#?'
        if not newresidue:
            newresidue='?'
        
        if not org_res:
            org_res='X'
        else:
            import DNAtool.mutation
            if DNAtool.mutation.three_to_one.has_key(org_res):
                org_res=DNAtool.mutation.three_to_one[org_res]
            else:
                org_res='?'

        import string
        if string.find(string.lower(parent),'select')!=-1:
            parent='P'
        if string.find(string.lower(newresidue),'select')!=-1:
            newresidue='X'
        else:
            import DNAtool.mutation
            newresidue=DNAtool.mutation.three_to_one[newresidue]
        import types
        if type(aa_number) is types.IntType:
            aa_number=int(aa_number)
        else:
            aa_number=str(int(str(aa_number[1:])))
        return '%s+%s:%s:%s:%s' %(parent,ChainID,aa_number,org_res,newresidue)


    def update_mut_selection(self, event=None):        
        """Update the data for the rest of the application"""
        
        if not hasattr(self,'structure_displayed'):
            self.structure_displayed=False
        
        newparent = self.parent.get()
        if newparent=='Select parent':
            newparent=None
        if newparent:
            aaseq=self.DB.get_protein_sequence(newparent)
            if aaseq:
                aaseq.sort()
                
                # Here we should change to allow for multiple chains                
                if newparent!=self.old_parent:
                    import string
                    number=string.split(aaseq[0][0],':')[1]
                    self.DT.resnum.set(number)
                    self.old_parent=newparent
                    
                    # If Yasara is started then delete the objects present                    
                    if self.yasara:
                        self.yasara.run('DelAll')
                        
                
                # Should we show the structure
                if DBActions.yasara != None:
                    pass
                #if int(self.show_struct.get())>0: # and self.yasara
                    
                    '''if not self.yasara_started:
                        print 'Setting up Yasara'
                        protein= self.DB[newparent]
                        DB = self.DB
                        if DB[newparent].has_key('Structure'):
                            import os
                            pdbname=os.path.join(os.getcwd(),'._pdbpipe')
                            fd=open(pdbname,'w')
                            fd.writelines(DB[protein]['Structure'])
                            fd.close()'''
             
                '''if self.yasara:
                    if self.yasara_started==1 and int(self.show_struct.get())>0:
                        
                        # Get Yasara to load the PDB file
                        
                        self.yasara.run('DelAll')
                        self.yasara.run('obj=LoadPDB %s' %pdbname)
                        self.yasara.run('Style Stick')
                        self.yasara.run('HideRes Hoh')
                        self.yasara.run('ColorAll Green')
                        self.yasara.run('HUD Off')
                        self.structure_displayed=True
                elif self.yasara and int(self.show_struct.get())>0 and self.structure_displayed:
                        #
                        # Rotate Yasara's view to the new residue and colour it red
                        #
                        number=self.DT.resnum.get()
                        self.yasara.run('ColorAll Green')
                        self.yasara.run('ColorRes %s,red' %(number))
                        if self.zoom_on_res.get()==1:
                            self.yasara.run('ZoomResidue %s,10' %str(number))'''
                        
                #self.update_mut_restext()
                name=self.suggest_name(self.parent.get(),
                                          self.DT.resnum.get(),
                                          self.DT.new_res.get(),
                                          self.DT.resname.get(),
                                          self.DT.chainID.get())
                self.new_name.set(name)
            else:
                import tkMessageBox
                tkMessageBox.showwarning('No protein sequence',
                                         'This parent has no protein sequence associated with it.\nI need a sequence to construct mutations!',
                                         parent=self.mut_win)
                return
        return
    
    def update_mut_restext(self,event=None):
        """When adding a mutant: Update the text in the orgres field"""
        
        # Find the original residue in the sequence        
        if not self.parent.get():
            return
     
        CID=self.DT.chainID.get()        
        # Update the to and from bounds for the scale        
        seq=self.DB[self.parent.get()]
        thischain=[]
        for resnum,aa in seq:
            if resnum.split(':')[0]==CID:
                thischain.append([resnum,aa])
        if thischain==[]:
            self.DT.resname.set('Unknown')
            return
        min_resnum=thischain[0][0].split(':')[1]
        self.DT.mutscale.configure(from_=min_resnum)
        max_resnum=thischain[-1][0].split(':')[1]
        self.DT.mutscale.configure(to=max_resnum)
        #
        #
        orgres=None
        aa_number=self.DT.resnum.get() # Get the residue number
        import string
        search_resnum='%s:%s' %(CID,string.zfill(aa_number,4))

        for resnum,aa in thischain:
            if resnum==search_resnum:
                orgres=aa
        if orgres:
            self.DT.resname.set(orgres)
        else:
            self.DT.resname.set('Unknown')
        return
        
