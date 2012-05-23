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

"""Routines for evaluating and entering primers"""

from Tkinter import *

class evaluate_primer:

    def __init__(self):
        """Dummy constructor"""
        import primer_database
        self.pDB=primer_database.primer_database(self,1)
        return

    def do_evaluate_primer(self,parent=None,i_parent=None,edit_primer_seq=None,
            edit_primer_descr=None, edit_primer_name=None):
        #
        # Import functions from primer_database
        #
        import primer_database
        self.pDB=primer_database.primer_database(i_parent=self,no_window=1)
        #
        # Who is the parent?
        #
        if not parent:
            parent=self.master
            self.i_parent=None
        else:
            self.i_parent=i_parent
        #
        # Open a window for entering primer data
        #
        self.eval_win=Toplevel()
        if edit_primer_seq:
            self.eval_win.title('Edit primer')
        else:
            self.eval_win.title('Add primer')
            
        self.eval_win.geometry('+%d+%d' %(parent.winfo_rootx()+50,
                                          parent.winfo_rooty()+50))
                                          
        
        # Create the string variables and attach them to entry widgets
        
        self.primername_var=StringVar()
        lab=Label(self.eval_win,text='Primer Name:')
        lab.grid(row=0,column=0,sticky='news')
        Ebox=Entry(self.eval_win,textvariable=self.primername_var,width=50)
        Ebox.grid(row=0,column=1,columnspan=2,sticky='news')
        
        self.eval_var=StringVar()
        lab=Label(self.eval_win,text='Sequence:')
        lab.grid(row=1,column=0,sticky='news')
        Ebox=Entry(self.eval_win,textvariable=self.eval_var,width=50)
        Ebox.grid(row=1,column=1,columnspan=2,sticky='news')
        
        self.descr_var=StringVar()
        desclab=Label(self.eval_win,text='Description:')
        desclab.grid(row=2,column=0,sticky='news')
        Ebox=Entry(self.eval_win,textvariable=self.descr_var,width=50)
        Ebox.grid(row=2,column=1,columnspan=2,sticky='news')
        
        #
        # Do recalculate and close buttons
        #
        row=3
        recalc=Button(self.eval_win,text='Recalculate',command=self.update_eval)
        recalc.grid(row=row,column=0,sticky='news')
        close=Button(self.eval_win,text='Close',command=self.close_eval)
        close.grid(row=row,column=1,sticky='news')
        #
        # Create button for adding to the primer library and set entry values if editing
        # a primer that is already present in the database
        #
        if edit_primer_seq:
            addp=Button(self.eval_win,text='Save primer',command=self.save_primer_from_edit)
            self.eval_var.set(edit_primer_seq)
            self.primername_var.set(edit_primer_name)
            self.descr_var.set(edit_primer_descr)
            self.edit_primer_name=edit_primer_name
        else:
            addp=Button(self.eval_win,text='Save primer',command=self.save_primer)
        addp.grid(row=row,column=2,sticky='news')
        
        #
        # Hairpin propensity
        row=row+1
        self.eval_hairpin=StringVar()
        self.eval_hairpin.set('No primer entered')
        l2=Label(self.eval_win,text='Hairpin propensity:')
        l2.grid(row=row,column=0,sticky='w')
        stat=Label(self.eval_win,textvariable=self.eval_hairpin)
        stat.grid(row=row,column=1,sticky='news')
        #
        # Primer-dimer
        #
        row=row+1
        self.eval_primerdimer=StringVar()
        self.eval_primerdimer.set('No primer entered')
        l2=Label(self.eval_win,text='Self complementarity:')
        l2.grid(row=row,column=0,sticky='w')
        stat=Label(self.eval_win,textvariable=self.eval_primerdimer)
        stat.grid(row=row,column=1,sticky='news')
        #
        # Tm
        #
        row=row+1
        self.eval_Tm=StringVar()
        self.eval_Tm.set('No primer entered')
        l2=Label(self.eval_win,text='Tm:')
        l2.grid(row=row,column=0,sticky='w')
        stat=Label(self.eval_win,textvariable=self.eval_Tm)
        stat.grid(row=row,column=1,sticky='news')
        #
        # Tm calculation method
        #
        self.eval_Tm_method=StringVar()
        self.eval_Tm_method.set(self.Tm_method)
        stat=Label(self.eval_win,textvariable=self.eval_Tm_method)
        stat.grid(row=row,column=2,sticky='news')
        #
        #
        #
        self.eval_win.bind('<KeyPress>',self.update_eval)
        self.eval_win.grab_set()
        Ebox.focus_set()
        self.update_eval()
        return self.eval_win

    def get_characteristics(self,primer,Tm_method):
        """Get the characteristics (hairpin prop, selfcomplementarity) of a primer"""
        #
        # Extract the DNA sequence
        #
        primer_seq=primer['sequence']
        import mutation
        ok,primer_seq=mutation.check_DNA(primer_seq)
        if ok:
            #
            # Hairpin propensity
            #
            prop,numscore=mutation.hairpin_propensity(primer_seq)
            if numscore < 8:
                score = 'OK'
            elif numscore >= 8:
                score = 'BAD'
            text_hairpin='%s (%2d matches)' %(score,numscore)
            #
            # Self complementarity
            #
            complementarity_dict,maxselfscore = mutation.self_complementarity_propensity(primer_seq)
            if maxselfscore < 8:
                text = 'OK'
            elif maxselfscore >= 8:
                text = 'BAD'
            text_selfcompl='%s (%2d matches)' %(text,maxselfscore)
            #
            # Melting temperature
            #
            Tm,mismatches,Tm_method_used=mutation.get_primer_Tm(DNAseq=primer['template_DNA'],primer=primer_seq,primer_start=primer['startpos'],method=Tm_method)
            return text_hairpin,text_selfcompl,Tm,Tm_method_used,mismatches
        return 'DNA sequenece not ok','DNA sequence not ok','DNA sequenece not ok','DNA sequence not ok'

    def update_eval(self,event=None):
        #
        # Calculate primer-dimer and secondary structure
        #
        import mutation
        primer_seq=self.eval_var.get()
        ok,primer_seq=mutation.check_DNA(primer_seq)
        #
        # Score the sequence
        #
        if ok:
            prop,numscore=mutation.hairpin_propensity(primer_seq)
            if numscore < 8:
                score = 'OK'
                # score = numscore
            elif numscore >= 8:
                score = 'BAD'
                #score = numscore
            self.eval_hairpin.set('%s (%2d matches)' %(score,numscore))
            #
            # Primer-Dimer
            #
            complementarity_dict,maxselfscore = mutation.self_complementarity_propensity (primer_seq)
            if maxselfscore < 8:
                text = 'OK'
            elif maxselfscore >= 8:
                text = 'BAD'
            self.eval_primerdimer.set('%s (%2d matches)' %(text,maxselfscore))
            #
            # Tm
            #
            # See if we can align the primer
            #
            Tm=None
            if self.data.has_key('DNAseq'):
                if self.data['DNAseq']:
                    import primer_alignment
                    A=primer_alignment.dummy_align_primer(self)
                    #
                    # Find all possible binding sites for the primer
                    #
                    sites=A.find_primer_binding_sites(primer_seq,self.data['DNAseq'])
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
                    Tm,mistmatches,Tm_method_used=mutation.get_primer_Tm(DNAseq=self.data['DNAseq'],
                                              primer=primer_seq,
                                              primer_start=best_position,
                                              method=self.Tm_method.get())
                    #
                    # Draw the primer on the screen in the best position
                    #
                    this_primer={}
                    this_primer['sequence']=primer_seq
                    this_primer['startpos']=best_position
                    #
                    # Draw the new primer
                    #
                    self.pDB.display_primer(this_primer)
            if not Tm:
                Tm,mistmatches,Tm_method_used=mutation.get_primer_Tm(DNAseq=None,
                                                         primer=primer_seq,
                                                         primer_start=None,
                                                         method=self.Tm_method.get())
            self.eval_Tm.set('%5.2f' %Tm)
            self.eval_Tm_method.set(Tm_method_used)
        else:
            self.eval_hairpin.set('Invalid DNA sequence')
            self.eval_primerdimer.set('Invalid DNA sequence')
            self.eval_Tm.set('Invalid DNA sequence')
        return

    def close_eval(self,event=None):
        """Close the primer evaluation window"""
        
        # Clear all graphics        
        if getattr(self,'pDB',None):
            self.pDB.clear_pDB_objects()
        
        # Close the window        
        self.eval_win.destroy()
        return

    def save_primer_from_edit(self, parent_window=None):
        """Save the primer in the DB when doing a primer Edit"""
        if not parent_window:
            parent_window=self.eval_win
        #
        # Check if the primer name has been changed first, and ask for confirmation of rename
        #        
        if self.primername_var.get() != self.edit_primer_name or not self.edit_primer_name:  
            import tkMessageBox
            ok = tkMessageBox.askyesno('Primer Name altered',
                                       'Rename primer and save?\n',
                                       parent=parent_window)
            if not ok:
                self.primername_var.set(self.edit_primer_name)
                return
            #rename primer based on value in entry widget    
            else:
                new_name=self.primername_var.get()
                if new_name:
                    #only rename if there isn't already a primer with the same name
                    if not self.data['primer_dict'].has_key(new_name):
                        self.data['primer_dict'][new_name]=self.data['primer_dict'][self.edit_primer_name].copy()
                        del self.data['primer_dict'][self.edit_primer_name]
                        self.edit_primer_name = new_name 
                        self.i_parent.new_name = new_name
                    #otherwise don't allow rename
                    else:
                        tkMessageBox.showwarning('Primer name already present',
                                     'Primer name already present.\nPlease choose another name',
                                     parent=parent_window)
                        return          

        #
        # Update the sequence and description fields for the primer
        #
        DNA=self.eval_var.get()
        # Validation check if primer is changed (maybe by mistake)
        if not DNA:            
            ok = tkMessageBox.askyesno('No sequence entered',
                                        'No sequence entered.\nDo you wish to save anyway?',
                                        parent=parent_window)
            if not ok:
                return
                
        import mutation
        ok,DNA_sequence=mutation.check_DNA(DNA)
        self.data['primer_dict'][self.edit_primer_name]['sequence']=DNA_sequence        
        self.data['primer_dict'][self.edit_primer_name]['description']=self.descr_var.get()
        # debug line
        #print 'Primer Info:',self.edit_primer_name,self.data['primer_dict'][self.edit_primer_name]
        #
        # Clear all graphics
        #
        if getattr(self,'pDB',None):
            self.pDB.clear_pDB_objects()
        #
        # Close the window
        #
        self.eval_win.destroy()
        return

    #
    # Saves primers to database
    #

    def save_primer(self,parent_window=None,openpdbWin=None,this_primer=None,calledby=None ):
        #
        # Get the parent window right
        #
        if not parent_window:
            parent_window=self.eval_win
        #
        # Save the primer to the primer library
        #
        if not self.data.has_key('primer_dict'):
            self.data['primer_dict']={}
        #
        # Find the next unique name
        #
        for x in range(1000):
            nextname='primer_%3d' %x
            if not self.data['primer_dict'].has_key(nextname):
                break
        #
        # Get a name for the primer
        # If calling from primer design just use little dialog, otherwise use
        #
        if calledby=='from design':
            print 'Saving primer from design window'
            done=None
            result=self.get_text_input(parent_window,[['name',20],['description',30]])
            if not result:
                #
                # User pressed Cancel
                #
                return None
            primer_name=result[0]
            description=result[1]
        else:
            primer_name=self.primername_var.get()
            description=self.descr_var.get()
            
        if not primer_name:
            import tkMessageBox
            tkMessageBox.showwarning('Invalid primer name',
                                     'I need a real primer name to store this primer',
                                     parent=parent_window)
            return
                
        print 'Setting name to',primer_name
        print 'Setting description to',description    
        
        # Store the primer        
        if not this_primer:
            DNA=self.eval_var.get()
            if not DNA:
                import tkMessageBox
                ok = tkMessageBox.askyesno('No sequence entered',
                                            'No sequence entered.\nDo you wish to save anyway?')
                if not ok:
                    return
            import mutation
            ok,sequence=mutation.check_DNA(DNA)
            if ok:
                self.data['primer_dict'][primer_name]={'description':description,'sequence':sequence,'startpos':None,
                                                       'hairpin_prop':self.eval_hairpin.get(),
                                                       'self-compl':self.eval_primerdimer.get(),
                                                       'introduced_sites':'Unknown'}
                self.primer_save_ok(primer_name,parent_window)
                #pass name of primer back to parent window for highlighting
                self.i_parent.new_name = primer_name
                parent_window.destroy()
                print 'Primer Info:',self.data['primer_dict'][primer_name]
            else:
                import tkMessageBox
                tkMessageBox.showwarning('Invalid primer','Cannot store an invalid primer',parent=parent_window)
                return
        else:
            
            #  If we already have all the data then save it            
            this_primer['description']=description
            self.data['primer_dict'][primer_name]=this_primer.copy()
            self.primer_save_ok(primer_name,parent_window)
        
        # Clear all the graphics        
        if getattr(self,'pDB',None):
            self.pDB.clear_pDB_objects()
            
        # If calling from mutagenic design, handle refresh of pdb here    
        if calledby=='from design' and openpdbWin!=None:
            openpdbWin.show_pDB_contents()
            openpdbWin.highlight_primer(primer_name)
            openpdbWin.display_details()

        return 

    def primer_save_ok(self,primer_name,parent_window):
        
        # Show info        
        import tkMessageBox
        tkMessageBox.showinfo('Primer saved','Primer %s saved to the Primer Database' %primer_name,parent=parent_window)
        
        # Mark that we need to save the project        
        self.data['Project saved']=None
        self.assess_status()
        return
        
