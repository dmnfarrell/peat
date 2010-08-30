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

# GUI for primer design

from Tkinter import *
#from tooltip import *
import Pmw

class primer_design_GUI:
    
    Tm_meth="Stratagene"
    def design_PCR_primer(self,parent=None):
        #
        # open window for designing pcr primers
        #
        self.pp_win=Toplevel()
        self.pp_win.title('PCR primers')
        self.pp_win.geometry('+100+500')
        Label(self.pp_win,text='Specify DNA fragment to amplify').grid(row=0,column=0)
        L1=Label(self.pp_win,text='First base:')
        L1.grid(row=1,column=0)
        self.first_base=IntVar()
        self.first_base.set(1)
        E1=Entry(self.pp_win,textvariable=self.first_base)
        E1.grid(row=1,column=1) 
        #
        #
        self.last_base=IntVar()
        self.last_base.set(1)
        L2=Label(self.pp_win,text='Last base:')
        L2.grid(row=2,column=0)
        E2=Entry(self.pp_win,textvariable=self.last_base)
        E2.grid(row=2,column=1)
        #
        # Length
        #
        self.length=IntVar()
        self.length.set(self.last_base.get()-self.first_base.get())
        Label(self.pp_win,text='Fragment length').grid(row=3,column=0)
        Label(self.pp_win,textvariable=self.length).grid(row=3,column=1)
        #
        #
        Label(self.pp_win,text='Specify fragment by click-and-drag in sequence window').grid(row=4,column=0,columnspan=3)
        Label(self.pp_win,text=' or by changing base numbers').grid(row=5,column=0,columnspan=3)
        #
        # Go button
        #
        go=Button(self.pp_win,text='Design primers',command=self.do_design_PCR)
        go.grid(row=6,column=0,columnspan=2,sticky='news')
        Button(self.pp_win,text='Cancel',command=self.cancel_design_PCR).grid(row=6,
                                                                              column=2,
                                                                              columnspan=1,
                                                                              sticky='news')
        #
        # Make the bindings
        #
        self.routine_call_DNAselection=self.update_PCR_primer
        return
    #
    # --------
    #

    def update_PCR_primer(self):
        """Update the field in the PCR primer design dialog"""
        self.first_base.set(self.data['DNA_selection']['start'])
        self.last_base.set(self.data['DNA_selection']['stop'])
        self.length.set(abs(self.last_base.get()-self.first_base.get()))
        return
    


    #
    # --------
    #

    def cancel_design_PCR(self,event=None):
        """Cancel the design of PCR primers"""
        self.pp_win.destroy()
        return

    #
    # ---------
    #

    def do_design_PCR(self):
        """Call the routines in mutation.py that design PCR primers"""
        return
        first_base=self.first_base.get()
        print 'FB',first_base
        import mutation
        self.results=mutation.PCR_primer(self.data['DNAseq'],65,first_base)
        self.listbox=Listbox(self.pp_win)
        self.listbox.grid(row=3,column=0)
        #
        self.text=Text(self.pp_win,width=20,bg='white')
        self.text.grid(row=3,column=3)
        #
        # Add items to the listbox
        #
        for primer in self.results.keys():
            self.listbox.insert(END,self.results[primer])
        #
        # Bind an action to the listbox
        #
        self.listbox.bind("<ButtonRelease-1>",self.display_PCR_primer)
        return

    #
    # ------
    #

    def display_PCR_primer(self,event=None):
        print self.listbox.curselection()
        selection=int(self.listbox.curselection()[0])
        primer_sel=self.results.keys()[selection]
        self.text.insert(END,primer_sel)
        return
        

    #
    # ----------------
    #

    def design_mutagenic_primer(self):
        """Design a mutagenic primer"""
        
        #
        # If the pDB is open, then make sure that we delete all the graphics
        #
        if getattr(self,'pDB_open',None):
            self.pDB_open.clear_pDB_objects()

            
        #
        # Open a new little window for specifying the PCR primers
        #
        self.primer_win=Toplevel()
        self.primer_win.transient(self.master)
        #self.primer_win.grab_set()
        self.primer_win.focus_set()
        self.primer_win.title('Mutagenic primer design')
        self.primer_win.geometry('+300+450')
        #self.master.wait_window(self.primer_win)

        H1=Label(self.primer_win,text='Specify desired characteristics of primer',
                 font='Times 16')
        H1.grid(row=0,column=0,columnspan=6,sticky='news')

        # ------------------------------------
        # Labels
        #
        row=2
        #
        # Specify the residue to mutate
        #
        row=self.specify_mut_res(window=self.primer_win,row=row)
        #
        # -------------------------
        # TM entry
        row=3
        l=Label(self.primer_win,text='Desired Tm')
        l.grid(row=row,column=0)
        self.TM_desired=DoubleVar()
        E1=Entry(self.primer_win,textvariable=self.TM_desired,width=6)
        E1.grid(row=row,column=1,sticky='news')
        self.TM_desired.set(65.0)
        
        # ------------------------------------
        # combo box for Tm methods
        #
        row=row+1
        def choseEntry(entry):
            self.Tm_meth=entry
        column=1
        methods = ("Simple", "Stratagene", "Basic") 
        Tmcombobox = Pmw.ComboBox(self.primer_win, label_text='Tm Method:', labelpos='wn',
                        listbox_width=24, dropdown=1,
                        selectioncommand=choseEntry,
                        scrolledlist_items=methods)
        #Tmtip = Pmw.Balloon(self)
        #Tmtip.bind(Tmcombobox, "Select TmMethod")
        Tmcombobox.grid(row=row,column=column,sticky='news')
        #Tmcombobox.pack(fill=BOTH, expand=1, padx=8, pady=8)
        Tmcombobox.selectitem(methods[1]) 
        
        
        # ----------------------------------------
        # Search for Tm
        #
        column=2
        l2=Label(self.primer_win,text='Design/remove restriction site')
        l2.grid(row=row,column=column+1,columnspan=1,sticky='nws')
        self.silent_mut=IntVar()
        E4=Checkbutton(self.primer_win,onvalue=1,offvalue=0,variable=self.silent_mut)
        E4.grid(row=row,column=column,sticky='news')
        self.silent_mut.set(1)         
        
        # ------------------------------------
        # Action buttons
        #
        row=row+1
        btn1=Button(self.primer_win,text='Cancel',command=self.close_primer_win)
        btn1.grid(row=row,column=0)
        
        btn2=Button(self.primer_win,text='Design primer',command=self._design_primer)
        #btn2tip = ToolTip(btn2, follow_mouse=1, text="Start search for primers")
        btn2.grid(row=row,column=1)
        #
        # Bindings
        #
        self.primer_win.bind('<Return>',self.update_restext)
        self.primer_win.bind('<Key>',self.update_restext)
        #
        # Create listbox for the results
        #
        row=row+1
        yscrollbar=Scrollbar(self.primer_win,orient='vertical',width=10)
        yscrollbar.grid(row=row,column=2,sticky='nws')
        #
        self.primer_result=Listbox(self.primer_win,bg='white',
                                   fg='black',
                                   height=10,width=30,yscrollcommand= yscrollbar.set)
        self.primer_result.grid(row=row,column=0,columnspan=2,sticky='news')
        yscrollbar.config(command=self.primer_result.yview)
        #
        # Detailed results window
        #
        self.detailed_results=Text(self.primer_win,background='white',
                                   foreground='black',state=NORMAL,exportselection=1,width=40,height=10)
        self.detailed_results.grid(row=row,column=3,columnspan=3,sticky='NWS')
        self.detailed_results.config(state=NORMAL)
        #
        # Action buttons for what to do
        #
        row=row+1
        column=0
        buttons=[['Save primer',self.save_primer_indesign],
                 ['Use primer (Save + send to EAT_DB)',self.use_primer],
                 ['Close',self.close_primer_win]]
        for button,action in buttons:
            x=Button(self.primer_win,text=button,command=action)
            x.grid(row=row,column=column)
            column=column+1
        #
        # Grab focus
        #
        self.primer_win.grab_set()
        E1.focus_set()

        #
        # Init view
        #
        self.update_restext()
        return

    #
    # ------
    #

    def save_primer_indesign(self):
        #
        # Save the selected primer to the primer database
        #
        if self.primer_result.curselection():
            selection=int(str(self.primer_result.curselection()[0]))
            name_selected=self.primer_order[selection]
            this_primer=self.foundprimers[name_selected]
            if getattr(self,'pDB_open',None):
                #pass a reference to the current primer db window	
                primerdbwin=self.pDB_open
                self.save_primer(self.primer_win,primerdbwin,this_primer,'from design')
            else:
                #If the pDB was not found open the window and lower it    
                import primer_database
                self.pDB_open=primer_database.primer_database(i_parent=self)
                self.pDB_open.pDB_win.lower()
                primerdbwin=self.pDB_open
                self.save_primer(self.primer_win,primerdbwin,this_primer,'from design')
            
        else:
            import tkMessageBox
            tkMessageBox.showwarning("No Selection",'Please create or select your primer',
                                     parent=self.primer_win)
            return
    def use_primer(self):
        
        return


    def specify_mut_res(self,window=None,row=None,command=None,
                        minres=None,maxres=None,ChainIDs=['A']):

        frame = Frame(window)
        frame.grid(row=row,column=0,columnspan=2)
        
        # Print the widgets for specifying the residue to be mutated        
        l2=Label(frame,text='Specify mutation:')
        l2.grid(row=row,column=0)

        # --------------------------------------------------------
        # New residue type
        # Make a long list of all the AA types you can mutate to
        
        self.new_res=StringVar()
        self.new_res.set('Select')
        self.newres_button=Menubutton(frame,textvariable=self.new_res,relief=RAISED)
        self.newres_menu=Menu(self.newres_button,tearoff=0)
        self.newres_button['menu']=self.newres_menu
        
        # New residue        
        column=1
        import mutation
        aas=mutation.three_to_one.keys()
        aas.sort()
        for aa in aas:
            self.newres_menu.add_radiobutton(label=aa,
                                             variable=self.new_res,
                                             value=aa,
                                             indicatoron=0,
                                             command=command)
        self.newres_button.grid(row=row,column=column,sticky='news')
        # Label
        lc=Label(frame,text='New AA')
        lc.grid(row=row-1,column=column)

        # Original Residue name        
        row=2
        column=2
        self.resname=StringVar()
        #E2=Entry(frame,textvariable=self.resname,width=4)
        E2=Label(frame,textvariable=self.resname)
        E2.grid(row=row,column=column,sticky='news')

        la=Label(frame,text='Org AA')
        la.grid(row=row-1,column=3)

        row=4
        # Chain Identifier
        column=0
        la=Label(frame,text='ChainID')
        la.grid(row=row-1,column=column)
        self.chainID=StringVar()
        self.chainID.set(ChainIDs[0])
        self.chainID_button=Menubutton(frame,textvariable=self.chainID,relief=RAISED)
        self.chainID_menu=Menu(self.chainID_button,tearoff=0)
        self.chainID_button['menu']=self.chainID_menu
        
        # add the chain ids        
        for CID in ChainIDs:
            self.chainID_menu.add_radiobutton(label=CID,
                                              variable=self.chainID,
                                              value=CID,
                                              indicatoron=0,
                                              command=command)
        self.chainID_button.grid(row=row,column=column,sticky='news')

        # Residue number
        # Get the start and end AA numbers
        
        column=1
        if not minres:
            try:
                min_resnum=self.data['ORF_selected']['aastart_number']
            except:
                min_resnum=1
        else:
            min_resnum=minres
        
        # Get the maximum residue number        
        if not maxres:
            max_resnum=min_resnum+len(self.data['ORF_selected']['aaseq'])-1
        else:
            max_resnum=maxres
        
        # Entry field for resnum        
        E3=Entry(frame,textvariable=self.resnum,width=10,justify='center')
        E3.grid(row=row,column=column,sticky='ns')
        # Label
        lb=Label(frame,text='Residue number',justify='center')
        lb.grid(row=row-1,column=column,columnspan=1,sticky='ns')
        #
        self.mutscale=Scale(frame,resolution=1,variable=self.resnum,orient='horizontal',
                            showvalue=0,command=self.update_restext,
                            from_=min_resnum,to=max_resnum,)
        self.mutscale.grid(row=row,column=column+1,sticky='news')
        self.resnum.set(min_resnum)        
        
        return row

    #
    # --------------
    #

    def update_restext(self,junk=None):
        #
        # Update the stringvar showing the old AA number
        #
        #frame=self.data['frame_sel']
        try:
            resnumber=self.resnum.get()-self.data['ORF_selected']['aastart_number']
            self.resname.set(self.data['ORF_selected']['aaseq'][resnumber])
            if resnumber>len(self.data['ORF_selected']['aaseq'])-1 or resnumber<0:
                # Out of range
                return
        except:
            #
            # User entered a non-number in the entry field
            #
            return
        #
        # Update the view of the sequence window
        #
        x,y=self.get_aa_pos_on_screen(resnumber+self.data['ORF_selected']['start'],0)
        #
        # We need to center on the aa preferably
        #
        center=max(0.0,x-(self.x_size-self.canvas_border_x)/2.0)
        #
        # Moveto works in fractions of the screen, so get the fraction
        #
        frac=center/self.canvas_x
        self.seqframe.xview('moveto', frac)
        #
        # Update the sequence window
        #
        self.update_sequence_window()
        return
    #
    # --------------
    #

    def close_primer_win(self):
        """Close the primer design window"""
        #
        # Delete all the graphics
        #
        if getattr(self,'pDB',None):
            self.pDB.clear_pDB_objects()
        #
        # Set the cancel flag
        #
        self.cancelled=1
        #
        # Close the window
        #
        self.primer_win.destroy()
        return
    #
    # -----------
    #
    
    def _design_primer(self):
        #
        # Clear text boxes
        #
        self.detailed_results.delete(1.0,END)
        self.primer_result.delete(0, END)
        #
        # Get an instance of the primer DB
        #
        if not getattr(self,'pDB',None):
            import primer_database
            self.pDB=primer_database.primer_database(self,1)
        #
        # Clear any previous details shown in the sequence window
        #
        self.pDB.clear_pDB_objects()
        #self.pDB.display_primer([],delete=1,only_delete=1)
        #
        # Get the melting temp
        #
        try:
            Tm_des=self.TM_desired.get()
            if Tm_des<40 or Tm_des>85:
                import tkMessageBox
                tkMessageBox.showwarning("Invalid Tm",'Please give a temperature in Celcius from 40 - 85',
                                         parent=self.primer_win)
                return
        except:
            import tkMessageBox
            tkMessageBox.showwarning("Invalid Tm",'Please give a temperature in Celcius from 40 - 85',
                                     parent=self.primer_win)
            return
        #
        # Get start and stop
        #
        try:
            AA_number=self.resnum.get()-self.data['ORF_selected']['aastart_number']
        except:
            #
            # Not a number
            #
            import tkMessageBox
            tkMessageBox.showwarning("Invalid AA number",'Please give a number from %d to %d'
                                     %(self.data['ORF_selected']['aastart_number'],
                                       self.data['ORF_selected']['aastart_number']+self.data['ORF_selected']['length']-1),
                                     parent=self.primer_win)
            return
        #
        # Check for out of range
        #
        if AA_number>self.data['ORF_selected']['length']-1 or AA_number<0:
            import tkMessageBox
            print self.data['ORF_selected']
            tkMessageBox.showwarning("AA number out of range",'Please give a number from %d to %d'
                                     %(self.data['ORF_selected']['aastart_number'],
                                       self.data['ORF_selected']['aastart_number']+self.data['ORF_selected']['length']-1),
                                       parent=self.primer_win)
                                       
            return
        #
        # Modify AA_number to get correct number in DNA sequence
        #
        AA_number_inORF=AA_number
        AA_number=AA_number+self.data['ORF_selected']['start']
        #
        new_AA=self.new_res.get()
        #
        # check that the user selected a new AA
        #
        if new_AA=='Select':
            import tkMessageBox
            tkMessageBox.showwarning("No new AA type selected",'You have to select a new residue type',
                                      parent=self.primer_win)
            return
        #
        # Warn if the new and old AA type are the same
        #
        old_AA=self.data['ORF_selected']['aaseq3'][AA_number_inORF]
        if new_AA==old_AA:
            import tkMessageBox
            if not tkMessageBox.askyesno("No AA change",
                                         "This mutation will not make a change in the amino acid sequence.\nContinue anyway?",
                                         parent=self.primer_win):
                return
        #
        # Should we look for a restriction site?
        #
        if self.silent_mut.get()==1:
            find_restriction_site=1
        else:
            find_restriction_site=None
        #
        # Adjust the DNA sequence according to frame
        #
        offset=self.data['ORF_selected']['frame']-1
        DNA_seq=self.data['DNAseq'][offset:]
        #
        # All ok
        #
        # Open a progress window
        #
        self.progress_win=Toplevel()
        self.progress_win.transient(self.primer_win)
        self.progress_win.focus_set()
        self.progress_win.title('Primer design progress')
        self.progress_win.geometry('+300+150')
        self.prog_x=300
        self.prog_can=Canvas(self.progress_win,bg='white',width=self.prog_x,height=180)
        self.prog_can.grid(row=0,column=0)
        #
        self.progress=IntVar()
        self.xl=Label(self.progress_win,textvariable=self.progress)
        self.xl.grid(row=1,column=0)
        self.box=None
        #
        self.master.update()
        #
        # Start the generation of primers
        #
        self.cancelled=None
        print 'AA_number',AA_number
        import mutation

        primers_results_dict, new_enzymes, enzymes_that_already_cut=mutation.exhaustive_research(DNA_seq,AA_number,
                                                                                                 new_AA,Tm_des,self.Tm_meth,
                                                                                                 enzyme_list=self.data['used_enzymes'],
                                                                                                 parent=self)
        if not primers_results_dict:
            self.progress_win.destroy()
            return
        #
        # Insert primers in the listbox
        #
        self.prog_can.create_text(5,105,text="Removing useless primers",anchor='nw')
        self.master.update()
        self.foundprimers={}
        count=0
        for primer_seq in primers_results_dict.keys():
            #
            # Filter the primers without new restriction sites if so specified
            #
            count=count+1
            name='Primer %4d' %count
            if not new_enzymes.has_key(primer_seq) and find_restriction_site:
                continue
            elif not new_enzymes.has_key(primer_seq):
                self.foundprimers[name]={}
                self.foundprimers[name]['introduced_sites']=[]
            else:
                self.foundprimers[name]={}
                self.foundprimers[name]['introduced_sites']=new_enzymes[primer_seq]
            #
            # Add the rest
            #
            self.foundprimers[name]['sequence']=primer_seq
            self.foundprimers[name]['Tm']=primers_results_dict[primer_seq]['Tm']
            #
            # Copy the other properties
            #
            for property in primers_results_dict[primer_seq].keys():
                self.foundprimers[name][property]=primers_results_dict[primer_seq][property]
        #
        # Delete primers that are not bringing any benefit
        #
        to_delete={}
        count=0
        for primer in self.foundprimers.keys():
            frac_done=float(count)/float(len(self.foundprimers.keys()))
            self.update_primer_progress(frac_done,4)
            count=count+1
            cuts=self.foundprimers[primer]['introduced_sites']
            mism1=self.foundprimers[primer]['mismatches']
            for primer2 in self.foundprimers.keys():
                if primer==primer2:
                    continue
                #
                # Compare with Other primer
                #
                if cuts==self.foundprimers[primer2]['introduced_sites']:
                    mism2=self.foundprimers[primer2]['mismatches']
                    if mism2<mism1:
                        to_delete[primer]=1
        #
        # Delete all primer that have been tagged
        #
        for delete in to_delete.keys():
            del self.foundprimers[delete]
        #
        # Fill the listbox, sort by number of mismatches
        #
        self.primer_result.delete(0, END)
        self.primer_names=self.foundprimers.keys()
        self.primer_names.sort()
        self.primer_order=[]
        for mismatch in range(20):
            for primer_name in self.primer_names:
                if self.foundprimers[primer_name]['mismatches']==mismatch:
                    self.primer_result.insert(END,primer_name+': %2d mism, Tm: %4.1f' %(mismatch,
                                                                                        self.foundprimers[primer_name]['Tm']))
                    self.primer_order.append(primer_name)
        #
        # Close the progress window
        #
        self.progress_win.destroy()
        self.master.update()
        #
        # Did we find any primers?
        #
        if len(self.primer_order)==0:
            import tkMessageBox
            tkMessageBox.showwarning('No primers found','I could not find any primers that match your criteria.\nAlter the query and try again.')
            return

        #
        # Bind button-1 presses to action
        #
        self.primer_result.bind("<Double-Button-1>",self.display_detailed_results)
        self.primer_result.bind("<ButtonRelease-1>",self.display_detailed_results)

        return

    #
    # ------------
    #

    def update_primer_progress(self,fraction,level=1):
        #
        #
        #
        y_pos=(level-1)*40+20
        bar_length=(self.prog_x-20)*fraction
        old_obj=self.box
        self.box=self.prog_can.create_polygon(10,y_pos,
                                              10+bar_length,y_pos,
                                              10+bar_length,y_pos+20,
                                              10,y_pos+20,fill='blue',
                                              outline='black')
        #
        # Text
        #
        txt='%4.1fx complete' %(100*fraction)
        import string
        txt=string.replace(txt,'x','%')
        self.progress.set(txt)
        if old_obj:
            self.prog_can.delete(old_obj)
        else:
            if level==1:
                bar_length=self.prog_x-20
                self.prog_can.create_polygon(10,y_pos,
                                             10+bar_length,y_pos,
                                             10+bar_length,y_pos+20,
                                             10,y_pos+20,fill='white',
                                             outline='black')
        self.master.update()
        return

    #
    # ------------
    #

    def display_detailed_results(self,junk=None):
        #
        # Figure out the selected primer and display in textbox
        #
        selection=int(str(self.primer_result.curselection()[0]))
        name_selected=self.primer_order[selection]
        this_primer=self.foundprimers[name_selected]
        self.detailed_results.delete(1.0,END)
        #
        # Insert the template DNA sequence
        #
        this_primer['template_DNA']=self.data['DNAseq']
        #
        # Get the alignment
        #
        this_primer['startpos']=self.pDB.align_primer(this_primer)
        #
        # Make sure that we have the updated the characteristics of the primer
        #
        import evaluate_primer
        EVAL=evaluate_primer.evaluate_primer()
        hairpin,selfcompl,Tm_inpos,Tm_method_used,mismatches=EVAL.get_characteristics(this_primer,self.Tm_method.get())
        this_primer['hairpin_prop']=hairpin
        this_primer['self-compl']=selfcompl
        #
        # Show all the info on the primer
        #
        self.detailed_results.insert(END,'Primer: %s \n' %name_selected)
        self.detailed_results.insert(END,'Length: %3d bases\n' %(len(this_primer['sequence'])))
        self.detailed_results.insert(END,"Forward 5' %s 3'\n" %this_primer['sequence'])
        #
        # Show the reverse complementary sequence
        #
        import mutation
        self.detailed_results.insert(END,"Reverse complementary: 5' %s 3'\n" %mutation.get_reverse_complementary(this_primer['sequence']))
        #
        # Show the characteristics of the primer
        #
        text='-----------------------\nCharacteristics\nTm (aligned): %5.2f (%s)\n' %(Tm_inpos,Tm_method_used)
        self.detailed_results.insert(END,text)
        text='Hairpin: %s, \nself-sim: %s, \nrestr. site differences: ' %(this_primer['hairpin_prop'],
                                                                        this_primer['self-compl'])
        self.detailed_results.insert(END,text)
        #
        # Show the recognition sequence for the enzyme(s)
        #
        unique_added,unique_removed,non_unique_added,non_unique_removed=self.pDB.get_primer_restriction_differences(this_primer)
        enz_specs=self.RS.enzymes_regexs
        inserted=None
        for enz in unique_added.keys():
            if enz_specs.has_key(enz):
                self.detailed_results.insert(END,enz+'(+)*, ')
                inserted=1
        #
        for enz in unique_removed.keys():
            if enz_specs.has_key(enz):
                self.detailed_results.insert(END,enz+'(-)*, ')
                inserted=1
        #
        for enz in non_unique_added.keys():
            if enz_specs.has_key(enz):
                self.detailed_results.insert(END,enz+'(+), ')
                inserted=1
        #
        for enz in non_unique_removed.keys():
            if enz_specs.has_key(enz):
                self.detailed_results.insert(END,enz+'(-), ')
                inserted=1
        #
        # If there were no differences in the restriction map, then write that
        #
        if inserted==None:
            self.detailed_results.insert(END,'None')
        #
        # Delete all graphic objects from last round
        #
        for obj in self.detailed_objs.keys():
            self.seqframe.delete(obj)
        #
        # Align the primer
        #
        this_primer['startpos']=self.pDB.align_primer(this_primer)
        #
        # Draw the new primer
        #
        self.pDB.display_primer(this_primer)
        return


    
    #
    # --------------
    #

    def about_primer_design(self):
        self.ab_win=Toplevel()
        self.ab_win.geometry('450x500+100+350')
        self.ab_win.title('About Primer design')
        text=['E A T D B','Enzyme Analysis Tool & DataBase',' ','A tool for analysing the effect of point mutations','on the catalytic and structural characteristics of proteins and enzymes','','Author: Jens Erik Nielsen','(C) Copyright 2005-2006 University College Dublin']
        row=0
        for line in text:
            tmp=Label(self.ab_win,text=line)
            tmp.grid(row=row,column=0)
            row=row+1
        return

if __name__=="__main__":
    print
    print 'Do not execute this file. Call DNAtool.py instead'
    print

