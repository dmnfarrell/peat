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

import sys,os,copy
import primer_alignment
from Tkinter import *
from PDBDumper import *
import Pmw

class primer_database(primer_alignment.align_primer,PDBDumper):

    """
    Primer database
    A single window that lists all primer names, - when clicked on, display sequence,
    Tm on its own, Tm when aligned to sequence, number of mismatches,
    changes to ORF, introduction of restriction site.
    Also display primer in the sequence window, aligned with  mismatches highlighted.
    """

    def __init__(self, i_parent, no_window=None, parentframe=None):
        self.parent=i_parent
        self.parentframe = parentframe
        #We just want to use the functions
        if no_window:
            return

        height=12
        details_width=50
        self.current_primers_shown=[]  #Track currently shown primer

        #if there is a parentframe given, we embed the primer db in it?
        '''if self.parentframe != None:
            print 'single window'
            self.pDB_win = Frame(self.parentframe,height=200)
            self.pDB_win.grid(row=6,column=0,rowspan=2,columnspan=4,padx=2,pady=2,sticky='news')
        #otherwise it's a separate window
        else:'''
        self.pDB_win=Toplevel()
        self.pDB_win.geometry('+200+480')
        self.pDB_win.title('Primer DataBase')
        self.pDB_win.protocol("WM_DELETE_WINDOW", self.close_pDB)

        #Pulldown menu
        self.balloon = Pmw.Balloon(self.pDB_win)
        menuBar = Pmw.MenuBar(self.pDB_win,
                            hull_relief = 'raised',
                            hull_borderwidth = 1,
                            hull_width=300,
                            balloon = self.balloon)
        menuBar.grid(row=0,column=0, columnspan=2, sticky='nw')

        menuBar.addmenu('File', 'primer database IO')
        menuBar.addmenuitem('File', 'command', 'Load primer database',
                command = self.load_primerDB,
                label = 'Load primer DB')
        menuBar.addmenuitem('File', 'command', 'Save primer database',
                command = self.save_primerDB,
                label = 'Save primer DB')
        menuBar.addmenuitem('File', 'command', 'Add primers from file',
                command = self.add_primers_to_DB,
                label = 'Add primers')
        menuBar.addmenuitem('File', 'command', 'Load from csv file',
                command = self.load_from_text,
                label = 'Load from csv file')
        menuBar.addmenuitem('File', 'command', 'Save to csv file',
                command = self.save_to_text,
                label = 'Save to csv file')
        menuBar.addmenuitem('File', 'command', 'Add primers from csv file',
                command = self.add_primers_from_text,
                label = 'Add primers from csv file')


        # Create listbox for the results
        row=1
        yscrollbar=Scrollbar(self.pDB_win,orient='vertical',width=14)
        yscrollbar.grid(row=row,column=2,sticky='nws')
        #
        self.primers=Listbox(self.pDB_win,bg='white',
                             fg='black',
                             height=height,width=30,yscrollcommand=yscrollbar.set,
                             selectmode=EXTENDED)
        self.primers.grid(row=row,column=0,columnspan=2,sticky='news',padx=2)
        yscrollbar.config(command=self.primers.yview)

        #Detailed results window
        detail_yscrollbar=Scrollbar(self.pDB_win,orient='vertical',width=12)
        detail_yscrollbar.grid(row=row,column=6,sticky='nws',padx=2)
        detail_xscrollbar=Scrollbar(self.pDB_win,orient='horizontal',width=12)
        detail_xscrollbar.grid(row=row+1,column=3,columnspan=3,sticky='news')
        self.details=Text(self.pDB_win,background='white',
                          foreground='black',
                          state=NORMAL,
                          exportselection=1,
                          width=details_width,
                          height=height,
                          yscrollcommand=detail_yscrollbar.set,
                          xscrollcommand=detail_xscrollbar.set,
                          wrap='word')
        detail_yscrollbar.config(command=self.details.yview)
        detail_xscrollbar.config(command=self.details.xview)
        self.details.grid(row=row,column=3,columnspan=3,sticky='NWS')
        self.details.config(state=DISABLED)
        self.pDB_win.grid_rowconfigure(row, weight=1)
        #self.pDB_win.grid_columnconfigure(3, weight=1)

        #use this variable to track the new primer name after add or edit
        self.new_name=''

        # Bottom panel consisting of Delete, Add new primer, Rename and close
        buttons=[['Delete',self.delete_primer],
                 ['Add',self.add_primer],
                 ['Edit',self.edit_primer],
                 ['Rename',self.rename_primer],
                 ['Apply primer',self.apply_primer],
                 ['Close',self.close_pDB]]
        row=row+2
        column=0
        for button,command in buttons:
            x=Button(self.pDB_win,text=button,command=command)
            x.grid(row=row,column=column,sticky='NEWS',padx=3,pady=2)
            column=column+1

        row=row+1
        self.find_type=StringVar()
        self.find_type.set('Find name')
        self.ftype_button=Menubutton(self.pDB_win,textvariable=self.find_type,relief=RAISED,width=16)
        self.ftype_menu=Menu(self.ftype_button,tearoff=0)
        self.ftype_button['menu']=self.ftype_menu

        # Other find types
        fts=['name','sequence','Tm','restriction site','mutation']
        for text in fts:
            text='Find '+text
            self.ftype_menu.add_radiobutton(label=text,
                                            variable=self.find_type,
                                            value=text,
                                            indicatoron=1)
        self.ftype_button.grid(row=row,column=0,sticky='nes',columnspan=2,padx=3,pady=2)

        self.find_var=StringVar()
        self.find_entry=Entry(self.pDB_win,background='white',textvariable=self.find_var,width=15)
        self.find_entry.grid(row=row,column=2,sticky='news',columnspan=2,padx=2,pady=2)
        Button(self.pDB_win,text='Find next',command=self.find_primer_inpDB).grid(row=row,column=4,sticky='news',columnspan=1,padx=2,pady=2)
        self.find_entry.bind('<KeyRelease>',self.find_primer_inpDB)
        self.find_solution=None

        self.orderbutton=Button(self.pDB_win,text='Create Order',command=self.write_current_primers)
        self.orderbutton.grid(row=row,column=5,sticky='news',columnspan=1,padx=2,pady=2)

        # Show the database contents
        self.show_pDB_contents()
        self.text_frame=None

        # Bind button-1 presses to action
        self.primers.bind("<Double-Button-1>",self.display_details)
        self.primers.bind("<ButtonRelease-1>",self.display_details)

        #get current path, if DB_Main is present then use project path
        if self.parent.parent:
            #if self.parent.parent.data['DBinstance']:
            if self.parent.parent.DB.data:
                self.path = os.path.dirname(self.parent.parent.DBfiledir)
                self.Dump = PDBDumper(self,self.path)
        #otherwise use current working directory of DNAtool
        else:
            self.path = os.getcwd()
            self.Dump = PDBDumper(self,self.path,noDB=1)

        # If main window is closed call close function
        #self.pDB_win.bind("<Destroy>", self.close_pDB )
        self.set_balloon_help()

        #recreate/clear primer order file for this session
        if os.path.isfile('current_order.txt'):
            os.unlink('current_order.txt')
        os.system('rm -rf current_order.txt')

        return


    def set_balloon_help(self):
        """Set help text for the popup balloons"""
        try:
            import Pmw
        except:
            import tkMessageBox
            tkMessageBox.showinfo("Missing component","To run this application you need to install Pmw (Python Megawidgets)\n or do 'Force update' to get the latest version of PEAT which includes Pmw'.")
            __EATDBversion__='0.000'

        self.balloon=Pmw.Balloon(self.pDB_win)
        help=[[self.orderbutton,'Write currently selected primers to order file']]
        for btn,txt in help:
            self.balloon.bind(btn,txt)
        return


    def find_primer_inpDB(self,event=None):
        """Find the primers in pDB that fulfills the criteria"""

        # Find all matches
        searchtext=(self.find_var.get()).upper()
        searchtype=self.find_type.get()

        solutions=[]
        for primer in self.primer_order:
            if searchtype=='Find name':
                if (primer.upper()).find(searchtext)!=-1:
                    solutions.append(primer)
            elif searchtype=='Find sequence':
                if (self.parent.data['primer_dict'][primer]['sequence']).upper().find(searchtext)!=-1:
                    solutions.append(primer)
            else:
                import tkMessageBox
                tkMessageBox.showwarning('not implemented yet','I never thought anyone would use this, so send me an email with a good reason why you need this, and I will fix it.\nJens.Nielsen@ucd.ie')
                return

        if not self.find_solution is None:
            self.find_solution=self.find_solution+1
            if self.find_solution>=len(solutions):
                self.find_solution=0
        else:
            self.find_solution=0

        # Highlight the solution and scroll the listbox

        self.primers.selection_clear(first=0,last=len(self.primer_order)-1)
        if len(solutions)>0 and searchtext!='':
            pname=solutions[self.find_solution]
            index=0
            for primer in self.primer_order:
                if primer==pname:
                    self.primers.see(index)
                if primer in solutions:
                    self.primers.selection_set(index)
                index=index+1
        return


    # Highlight a specified single primer

    def highlight_primer(self,pname):
        """Highlight and update details for the given primer after addition/editing"""
        self.primers.selection_clear(first=0,last=1)
        if pname!='':
            index=0
            for primer in self.primer_order:
                if primer==pname:
                    self.primers.see(index)
                    self.primers.selection_set(index)
                    return
                else:
                    index=index+1

    def primer_sort(self,primer1,primer2):
        """This function is for sorting primers"""
        p1=primer1.lower()
        p2=primer2.lower()
        if p1<p2:
            return -1
        elif p1>p2:
            return 1
        return 0

    def show_pDB_contents(self):

        # Show the database, sort by name
        self.primers.delete(0, END)
        if self.parent.data.has_key('primer_dict'):
            self.primer_names=self.parent.data['primer_dict'].keys()
        else:
            self.primer_names=[]
        self.primer_names.sort(cmp=self.primer_sort)
        self.primer_order=[]
        for primer in self.primer_names:
            self.primers.insert(END,primer)
            self.primer_order.append(primer)
        return

    def display_details(self,event=None):
        """Show all details on a selected primer. If we have a parent DNA sequence then
        show the alignment to that one"""

        # Figure out the selected primer and display in textbox
        tmp_selection=self.primers.curselection()
        selection=[]
        for num in tmp_selection:
            selection.append(int(num))

        # Display all of the primers
        primer_sites=[]
        first=1
        for primer_num in selection:
            site=self.display_single_primer(primer_num,focus=first,delete=first)
            first=None
            primer_sites.append(site)

        # Select the maximum length DNA sequence spanned by the primers
        if primer_sites!=[] and primer_sites != None:
            minstart=999
            maxend=-9999
            for start,end in primer_sites:
                minstart=min(start,minstart)
                maxend=max(end,maxend)

            # Highlight the sequence
            self.parent.clear_selection()
            self.parent.data['DNA_selection']={}
            self.parent.data['sel_objs']={}
            self.parent.data['DNA_selection']['start']=minstart
            self.parent.data['DNA_selection']['stop']=maxend
            self.parent.mark_base_selection(self.parent.data['DNA_selection']['start'],self.parent.data['DNA_selection']['stop'])
            self.parent.DNAfragment_counter.set(abs(self.parent.data['DNA_selection']['stop']-self.parent.data['DNA_selection']['start']))
            if getattr(self.parent,'routine_call_DNAselection',None):
                (self.parent.routine_call_DNAselection)()
        return

    def display_single_primer(self,selection,focus,delete):
        """Display a single primer - focus on it if focus is true, and delete
        previous graphics and text objects if delete is true"""
        name_selected=self.primer_order[selection]
        #print 'SELECTION=',selection
        this_primer=self.parent.data['primer_dict'][name_selected]
        self.current_primers_shown.append(this_primer)
        self.parent.primer_displayed=1  #tells DNAtool that a primer is currently shown
        self.details.config(state=NORMAL)
        if delete:
            self.details.delete(1.0,END)

        # Align the primer
        # If we do not have a DNA sequence then we cannot do it

        this_primer['startpos']=None
        if self.parent.data.has_key('DNAseq'):
            if self.parent.data['DNAseq']:

                # Find all possible binding sites for the primer
                sites=self.find_primer_binding_sites(this_primer['sequence'],self.parent.data['DNAseq'])

                # Print the number of binding sites
                best_score=0
                best_position=None
                scores=[]
                first_neg=1
                for position,score in sites:
                    scores.append(score)

                    # Find the best position
                    if score>best_score:
                        best_score=score
                        best_position=position

                # Set the primer start for the primer characteristics
                this_primer['startpos']=best_position
                this_primer['template_DNA']=self.parent.data['DNAseq']

        # Show all the info on the primer
        self.details.tag_config('n', foreground='blue')
        self.details.insert(END,'Primer: %s \n' %name_selected, 'n')
        self.details.insert(END,'Description: %s \n' %(this_primer['description']))
        self.details.insert(END,'Length: %3d bases\n' %(len(this_primer['sequence'])))
        self.details.insert(END,"Forward 5' %s 3'\n\n" %this_primer['sequence'])

        # Show the reverse complementary sequence
        import mutation
        self.details.insert(END,"Reverse complementary: 5' %s 3'\n" %mutation.get_reverse_complementary(this_primer['sequence']))

        # Make sure that we have the updated the characteristics of the primer
        if not this_primer.has_key('template_DNA'):
            return
        import evaluate_primer
        EVAL=evaluate_primer.evaluate_primer()
        hairpin,selfcompl,Tm_inpos,Tm_method_used,mismatches=EVAL.get_characteristics(this_primer,self.parent.Tm_method.get())
        this_primer['hairpin_prop']=hairpin
        this_primer['self-compl']=selfcompl
        this_primer['introduced_sites']='Unknown'

        # Show the characteristics of the primer
        text='-----------------------\nCharacteristics\nTm (in aligned position): %5.2f (%s)\n' %(Tm_inpos,Tm_method_used)
        self.details.insert(END,text)
        text='Hairpin: %s, \nself-sim: %s, \nrestr. site differences: ' %(this_primer['hairpin_prop'],
                                                                        this_primer['self-compl'])
        self.details.insert(END,text)
        #
        # Show the recognition sequence for the enzyme(s)
        #
        self.details.insert(END,'Unique sites are marked with a "*"\n')
        unique_added,unique_removed,non_unique_added,non_unique_removed=self.get_primer_restriction_differences(this_primer)
        enz_specs=self.parent.RS.enzymes_regexs
        inserted=None
        for enz in unique_added.keys():
            if enz_specs.has_key(enz):
                self.details.insert(END,enz+'(+)*, ')
                inserted=1

        for enz in unique_removed.keys():
            if enz_specs.has_key(enz):
                self.details.insert(END,enz+'(-)*, ')
                inserted=1
        for enz in non_unique_added.keys():
            if enz_specs.has_key(enz):
                self.details.insert(END,enz+'(+), ')
                inserted=1

        for enz in non_unique_removed.keys():
            if enz_specs.has_key(enz):
                self.details.insert(END,enz+'(-), ')
                inserted=1

        # If there were no differences in the restriction map, then write that
        if inserted==None:
            self.details.insert(END,'None')

        # Delete all graphic objects from last round
        if delete:
            if not getattr(self.parent,'detailed_objs',None):
                self.parent.detailed_objs={}
            for obj in self.parent.detailed_objs.keys():
                self.parent.seqframe.delete(obj)
                if getattr(self.parent,'tempsites',None):
                    self.parent.tempsites.remove(obj)
        #
        # Print guidelines for this primer
        #
        scores.sort()
        diff=scores[-1]-scores[-2]
        self.details.insert(END,'\nDifference in # of matches between two best sites: %2d\n' %diff)
        if diff<5:
            self.details.insert(END,'WARNING: Primer is not unique!\n')
        if scores[-1]!=len(this_primer['sequence']):
            self.details.insert(END,'\nWARNING: No perfectly matching binding site.\n')

        # Print the number of mismatches as a control

        if this_primer['startpos']>0:
            self.details.insert(END,'\nDisplaying position: %4d on forward strand\n' %this_primer['startpos'])
        else:
            self.details.insert(END,'\nDisplaying position: %4d on reverse strand\n' %(-this_primer['startpos']))
        self.details.insert(END,'Number of mismatches: %d\n\n' %mismatches)

        # Display the primer

        match=self.display_primer(this_primer,focus)

        # Show the binding site

        self.details.insert(END,'\n==================================\n')
        site=sites[0][0]
        if site<0:
            site=len(this_primer['sequence'])+site
        #
        DNA_stretch=self.parent.data['DNAseq'][site:site+len(this_primer['sequence'])]
        self.details.insert(END,'Primer %s\n       %s\nDNASeq %s\n' %(this_primer['sequence'],match,DNA_stretch))
        #
        # Print summary of all binding sites
        #
        self.details.insert(END,'\n---------------------\nSummary of best binding sites\n')
        self.details.insert(END,'Forward strand:\n')
        first_neg=1
        for site,score in sites:
            if site<0:
                site=site*-1
                if first_neg:
                    first_neg=None
                    self.details.insert(END,'Reverse compl. strand:\n')
            self.details.insert(END,'Position: %4d, matches: %2d\n' %(site,score))
        #
        # All done
        #
        self.details.insert(END,'\n')
        self.details.config(state=DISABLED)
        #
        # Return the sites for use in highlighting
        #
        start_position=self.get_real_primer_position(this_primer)
        return [start_position,start_position+len(this_primer['sequence'])]

    #
    # Refresh display of the primer in the sequence window. Done when formatting
    # is changed in sequence window

    def refresh_primer(self):
        """Refresh display of the primer in the sequence window"""

        #get current selection and just refresh it with new font/size
        #selection=self.primers.curselection()
        #find if there are any primers shown in sequence window
        #c=self.parent.seqframe
        #currprimers = c.find_withtag('primer')
        #print 'CURRENT PRIMERS',currprimers
        #c.delete(currprimers)

        #redisplay the stored list of currently displayed primers
        for currprimer in self.current_primers_shown:
            self.display_primer(currprimer)
            #also refresh current base selection
            if self.parent.data.has_key('DNA_selection'):
                self.parent.mark_base_selection(self.parent.data['DNA_selection']['start'],
                                                self.parent.data['DNA_selection']['stop'])
        return
    #
    # -----------
    #

    def display_primer(self,this_primer,delete=1,only_delete=None,focus=1):
        """Display the primer in the sequence window"""

        font = self.parent.getCurrentFont()
        # Keep track of the objects
        if not getattr(self,'primer_objs',None):
            self.primer_objs={}
        if delete or only_delete:
            for obj in self.primer_objs.keys():
                self.parent.seqframe.delete(obj)
            self.primer_objs={}
            if only_delete:
                return

        # Do we have a primer?
        if this_primer['sequence']=='' or not this_primer['sequence']:
            return

        # Draw the new primer
        mismatches,match,objs=self.draw_primer(this_primer,lvl=self.parent.maxseqlevel)
        for obj in objs.keys():
            self.primer_objs[obj]=1

        # Get the restriction map differences
        unique_added,unique_removed,non_unique_added,non_unique_removed=self.get_primer_restriction_differences(this_primer)

        # Plot the thing
        self.parent.plot_restriction_sites(unique_added,colour='darkgreen',direction='down',
                                           add_text='(+)',temporary=1,delete_temporary=1)
        self.parent.plot_restriction_sites(unique_removed,colour='red',direction='down',
                                           add_text='(-)',temporary=1,delete_temporary=None)
        self.parent.plot_restriction_sites(non_unique_added,colour='darkgreen',direction='down',
                                           add_text='(+)',temporary=1,delete_temporary=None,underline_unique=None)
        self.parent.plot_restriction_sites(non_unique_removed,colour='darkred',direction='down',
                                           add_text='(-)',temporary=1,delete_temporary=None,underline_unique=None)

        # See if there are any AA changes and show (only if we have an ORF)

        self.mutations=[]
        if self.parent.data.has_key('ORF_selected'):
            # Find any differences in the sequence
            frame=self.parent.data['ORF_selected']['frame']-1
            wtseq=self.parent.data['ORF_selected']['aaseq3']
            new_DNA=self.apply_primer_to_DNA(this_primer,self.parent.data['DNAseq'],this_primer['startpos'])
            import mutation
            AA_seqs3,AA_seqs1=mutation.translate(new_DNA)

            # Get only the part of sequence that's defined in the ORF
            ORF_start=self.parent.data['ORF_selected']['start']
            AA_seqs3[frame]=AA_seqs3[frame][ORF_start-1:]
            for count in range(min(len(AA_seqs1[frame]),len(wtseq))):
                if AA_seqs3[frame][count]!=wtseq[count]:
                    newres=AA_seqs3[frame][count]
                    position=count+1
                    #print 'position',position
                    x,y_junk=self.parent.get_aa_pos_on_screen(position+ORF_start-2,frame+2)
                    colour='red'
                    y=(y_junk-25)*self.parent.y_scale

                    obj=self.parent.seqframe.create_text(x+2,y,
                                                         text=newres,
                                                         font=font,
                                                         anchor='w',fill=colour)
                    self.primer_objs[obj]=1
                    self.mutations.append(':'+wtseq[count]+':'+str(count+self.parent.data['ORF_selected']['aastart_number'])+':'+AA_seqs3[frame][count])

        # Update the view of the sequence window
        if focus:
            position=this_primer['startpos']
            if this_primer['startpos']<0:
                position=len(self.parent.data['DNAseq'])-(abs(position)+len(this_primer['sequence']))
            x,y=self.parent.get_base_pos_on_screen(position+len(this_primer['sequence'])/2)

            # We need to center on the aa preferably
            center=max(0.0,x-(self.parent.x_size-self.parent.canvas_border_x)/2.0)
            # Moveto works in fractions of the screen, so get the fraction
            frac=center/self.parent.canvas_x
            self.parent.seqframe.xview('moveto', frac)
            self.parent.seqframe.yview('moveto', 0.2)  #middle of scrollregion
        return match


    def get_primer_restriction_differences(self,this_primer):
        """Find the impact of this primer of the restriction map of the parent DNA"""

        # Do a restriction digest with the new DNA sequence and display any added/removed sites
        new_DNA=self.apply_primer_to_DNA(this_primer,self.parent.data['DNAseq'],this_primer['startpos'])
        #
        self.parent.select_enzymes()
        new_cut=self.parent.RS.get_restriction_sites(new_DNA,self.parent.data['used_enzymes'])

        # Wild type restriction digest
        self.parent.select_enzymes()
        wt_cut=self.parent.RS.get_restriction_sites(self.parent.data['DNAseq'],
                                                    self.parent.data['used_enzymes'])

        # Find the difference between the two sets of digests
        enzymes=new_cut.keys()
        #to_be_plotted={}
        unique_added={}
        non_unique_added={}
        non_unique_removed={}

        for enz in enzymes:
            if wt_cut.has_key(enz):
                if wt_cut[enz]!=new_cut[enz]:
                    # Check if we removed non-unique cuts
                    for cut in wt_cut[enz]:
                        if not cut in new_cut[enz]:
                            # There's a cut missing for this enzyme when the primer is applied
                            if not non_unique_removed.has_key(enz):
                                non_unique_removed[enz]=[]
                            non_unique_removed[enz].append(cut)

                    # Check if we added non-unique cuts
                    for cut in new_cut[enz]:
                        if not cut in wt_cut[enz]:
                            # We added this one
                            if not non_unique_added.has_key(enz):
                                non_unique_added[enz]=[]
                            non_unique_added[enz].append(cut)

            else:
                # Added completely new site
                unique_added[enz]=new_cut[enz][:]

        # Check the other way
        unique_removed={}
        for enz in wt_cut.keys():
            if not enz in enzymes:
                unique_removed[enz]=wt_cut[enz][:]
        return unique_added, unique_removed,non_unique_added,non_unique_removed


    def delete_primer(self):
        """Delete a primer from the database"""

        self.clear_pDB_objects()

        # Delete a primer
        selection=int(str(self.primers.curselection()[0]))
        name_selected=self.primer_order[selection]
        print 'in primer_dict, deleting',name_selected
        del self.parent.data['primer_dict'][name_selected]

        # Update the view
        self.details.delete(1.0,END)
        self.show_pDB_contents()
        self.parent.projectChanged()
        return

    def add_primer(self):
        """Open the evaluate primer window to add a primer manually"""

        self.clear_pDB_objects()
        # Add a primer
        win = self.parent.do_evaluate_primer(self.pDB_win,self)
        self.pDB_win.wait_window(win)

        # Update view
        self.show_pDB_contents()
        #print 'PRIMER_NAME=',self.new_name
        self.highlight_primer(self.new_name)
        self.display_details()
        return

    def edit_primer(self):

        """Edit a primer, using the current selection to locate the primer name"""

        try:
            selection=int(str(self.primers.curselection()[0]))
        except:
           import tkMessageBox
           tkMessageBox.showwarning('No primer selected',
                                     'Please select a primer in the listbox',
                                     parent=self.pDB_win)
           return

        selection=int(str(self.primers.curselection()[0]))
        name_selected=self.primer_order[selection]

        # Undisplay the original primer since do_evaluate_primer will show the primer while editing

        if not getattr(self,'primer_objs',None):
            self.primer_objs={}
        for obj in self.primer_objs.keys():
            self.parent.seqframe.delete(obj)
            del self.primer_objs[obj]

        # Call do_evaluate primer

        this_primer=self.parent.data['primer_dict'][name_selected]
        win=self.parent.do_evaluate_primer(self.pDB_win,self,
                                           edit_primer_seq=this_primer['sequence'],
                                           edit_primer_descr=this_primer['description'],
                                           edit_primer_name=name_selected)
        self.pDB_win.wait_window(win)
        #
        # Update view
        #
        self.show_pDB_contents()
        print 'PRIMER_NAME=',self.new_name
        self.highlight_primer(self.new_name)
        self.display_details()
        return


    def rename_primer(self):
        #
        # Rename a primer
        #
        try:
           selection=int(str(self.primers.curselection()[0]))
        except:
           import tkMessageBox
           tkMessageBox.showwarning('No primer selected',
                                     'Please select a primer in the listbox',
                                     parent=self.pDB_win)
           return

        name_selected=self.primer_order[selection]
        #
        # Get the new name
        #
        import tkSimpleDialog
        new_name=tkSimpleDialog.askstring('Rename primer','Enter new name',
                                          initialvalue=name_selected,
                                          parent=self.pDB_win)
        if new_name:
            if not self.parent.data['primer_dict'].has_key(new_name):
                self.parent.data['primer_dict'][new_name]=self.parent.data['primer_dict'][name_selected].copy()
                del self.parent.data['primer_dict'][name_selected]
        #
        # Update view
        #
        self.show_pDB_contents()
        self.highlight_primer(new_name)
        self.display_details()
        return

    #
    # --------------
    #

    def align_primer(self,this_primer):
        """Align the primer to the template DNA"""
        #
        # Find out where the primer is aligned
        #
        sites=self.find_primer_binding_sites(this_primer['sequence'],self.parent.data['DNAseq'])
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
        return best_position

    #
    # -------
    #

    def apply_primer(self):
        """Apply a primer to the DNA sequence"""
        selection=None
        #
        # We only accept a single primer
        #
        if len(self.primers.curselection())>1:
            import tkMessageBox
            tkMessageBox.showinfo('Apply primer',
                                  'You can only apply a single primer at the time',
                                  parent=self.pDB_win)
            return
        #
        # First a warning
        #
        import tkMessageBox
        if not tkMessageBox.askyesno('Apply primer',
                                     'I will apply the primer to the main DNA sequence\n\nDo you want to continue?.',
                                     parent=self.pDB_win):
            print 'I am not applying the primer!'
            return
        #
        # Apply the primer
        #
        #
        # First find the best binding site
        #
        selection=int(str(self.primers.curselection()[0]))
        name_selected=self.primer_order[selection]
        this_primer=self.parent.data['primer_dict'][name_selected]
        this_primer['startpos']=self.align_primer(this_primer)
        #
        # Apply the primer
        #
        new_DNA=self.apply_primer_to_DNA(this_primer,self.parent.data['DNAseq'],this_primer['startpos'])
        #
        # Clear all EAT_DB specific records
        #
        self.parent.data['DNAseq_status']='ALTERED'
        self.parent.data['DNAseq_mutations']=self.mutations[:]
        self.parent.data['DNAseq']=new_DNA
        #
        # Update the view
        #
        self.parent.update_sequence_window()
        self.display_details()
        return

    #
    # -----------
    #

    def apply_primer_to_DNA(self,this_primer,DNA_seq,position):
        #
        # Construct new DNA sequence
        #
        # If we have a negative position, then it's a reverse strand primer
        #
        if position<0:
            import mutation
            position=abs(int(position))
            repl_seq=mutation.get_reverse_complementary(this_primer['sequence'])
            position=len(DNA_seq)-(position+len(repl_seq))
        else:
            repl_seq=this_primer['sequence']
        #
        # Construct the new DNA
        #
        new_DNA=DNA_seq[:position]+repl_seq+DNA_seq[position+len(repl_seq):]
        return new_DNA

    #
    # ----------
    #

    def clear_pDB_objects(self):
        """Clear all objects that were displayed temporarily by primer_database functions
        Also clear the self.mutations variable"""
        #
        # Delete all graphics that we constructed
        #
        if not getattr(self.parent,'detailed_objs',None):
            self.parent.detailed_objs={}
        for obj in self.parent.detailed_objs.keys():
            self.parent.seqframe.delete(obj)
        self.parent.details_objs={}
        #
        # Remove all the temporary restriction sites
        #
        if getattr(self.parent,'temp_objs',None):
            for obj in self.parent.temp_objs.keys():
                self.parent.seqframe.delete(obj)
            self.parent.temp_objs={}
        #
        # Another dict
        #
        if getattr(self,'primer_objs',None):
            for obj in self.primer_objs.keys():
                self.parent.seqframe.delete(obj)
            self.primer_objs={}
        #also remove rects and lines for temp sites
        self.parent.seqframe.delete('templabelrect')
        self.parent.seqframe.delete('templine')
        #
        # Clear the mutations variable
        #
        self.mutations=[]
        return

    #
    # ----------
    #

    def close_pDB(self,event=None):
        """Close the primer database and clean up the display"""
        #self.Dump.doDump(event)
        print 'closing pdb window'
        self.clear_pDB_objects()
        self.parent.primer_displayed==0
        #
        # Close window
        #
        self.pDB_win.destroy()
        #
        # Reset the variable in the parent
        #
        self.parent.pDB_open=None
        return


    def load_primerDB(self,overwrite=1):
        """Load a new primer database"""
        import tkFileDialog, os
        initialdir=os.getcwd()
        filename=tkFileDialog.askopenfilename(defaultextension='primerDB',
                                              initialdir=initialdir,
                                              filetypes=[("Primer database files","*.primerDB"),
                                                         ("All files","*.*")])
        if not filename:
            return

        # Load a simple pickled file (for now)

        import pickle
        fd=open(filename,'rb')
        pdict=pickle.load(fd)
        fd.close()

        # Should we overwrite the existing primers?
        if overwrite:
            import tkMessageBox
            '''ans = tkMessageBox.askyesno("Overwrite primer database?",
                     "This will overwrite all primers you have in the primer datbase.\nAre you sure you want to continue?")'''
            from PEATDB.Dialogs import askyesnocancel
            a = askyesnocancel(title='Overwrite primer database?',
                               message='This will overwrite all primers you have in the primer datbase.\nAre you sure you want to continue?',
                               parent=self.pDB_win)
            ans = a.result
            print ans
            if ans == 'yes':
                self.parent.data['primer_dict']=pdict.copy()
                self.show_pDB_contents()
        return pdict


    def load_from_text(self,overwrite=1):
        """Load a new primer database from a csv text file"""
        import tkFileDialog, os
        initialdir=os.getcwd()
        filename=tkFileDialog.askopenfilename(defaultextension='pdb.csv',
                                              initialdir=initialdir,
                                              filetypes=[("Text files","*.csv"),
                                                         ("All files","*.*")])
        if not filename:
            return
        pdict={}
        import csv
        reader = csv.reader(open(filename, "rb"))
        for row in reader:
            try:
                primer=row[0]
                pdict[primer]={'description':row[1],'sequence':row[2]}
            except:
                print 'failed to load'
        #print pdict

        #
        # Should we overwrite the existing primers?
        #
        if overwrite:
            import tkMessageBox
            if tkMessageBox.askyesno("Overwrite primer database?",
                                     "This will overwrite all primers you have in the primer datbase.\nAre you sure you want to continue?"):
                self.parent.data['primer_dict']=pdict.copy()
                self.show_pDB_contents()
        return pdict

    #
    # ----
    #

    def save_primerDB(self):
        """Save the primers to a file"""
        import tkFileDialog, os
        savedir=os.getcwd()
        filename=tkFileDialog.asksaveasfilename(defaultextension='.primerDB',
                                                initialdir=savedir,
                                                filetypes=[("Primer database files","*.primerDB"),("All files","*.*")])
        if not filename:
            return
        #
        # Write the file
        #
        primers=self.parent.data['primer_dict'].copy()
        fd=open(filename,'wb')
        import pickle
        pickle.dump(primers,fd)
        fd.close()
        return


    #
    # ----
    #

    def save_to_text(self,filename=None,separator=',',nodesc=None):
        """Save the primers to a csv text file. Same format as dump file"""
        import tkFileDialog, os
        savedir=os.getcwd()
        if filename==None:
            filename=tkFileDialog.asksaveasfilename(defaultextension='.pdb.csv',
                                                    initialdir=savedir,
                                                    filetypes=[("Text files","*.csv"),("All files","*.*")])

        if not filename:
            return
        #
        # Write the file
        #
        primers=self.parent.data['primer_dict'].copy()
        fd=open(filename,'wb')

        #this code reused from the dumper class, some duplication
        if nodesc==None:
            HeaderList = ['description', 'sequence']
        else:
            HeaderList = ['sequence']
        #ActualHeaderList = ['name','description', 'sequence']
        DumpPDB = []
        #DumpPDB.append(ActualHeaderList)
        for entry in primers:
            tmp = []
            tmp.append(entry)
            for head in HeaderList:
                tmp.append(primers[entry][head])
            DumpPDB.append(tmp)
        result = DumpPDB

        DFILE = filename
        try:
            fd=open(DFILE,'w')
            for line in result:
                str = ''
                for S in line:
                    str = str+S+separator
                fd.write(str+'\n')
            fd.close()
            print "File written ",DFILE
        except:
            print "error: could not write file ",DFILE

        fd.close()
        return

    #
    # ----
    #

    def add_primers_to_DB(self):
        """Add primers from a file"""
        pdict=self.load_primerDB(overwrite=None)
        primers=pdict.keys()
        primers.sort()
        import copy
        for primer in primers:
            if self.parent.data['primer_dict'].has_key(primer):
                import tkMessageBox, copy
                if tkMessageBox.askyesno("Overwrite primer?",
                                         "The file contains a primer named '%s', but your primer database already\ncontains a primer with this name.\nDo you want to overwrite your existing primer with the new one?" %primer):
                    self.parent.data['primer_dict'][primer]=copy.deepcopy(pdict[primer])
            else:
                self.parent.data['primer_dict'][primer]=copy.deepcopy(pdict[primer])
        self.show_pDB_contents()
        return


    #
    # Add primers from a csv file
    #
    def add_primers_from_text(self):
        """Add primers from a csv file"""
        pdict=self.load_from_text(overwrite=None)
        primers=pdict.keys()
        primers.sort()
        import copy
        for primer in primers:
            if self.parent.data['primer_dict'].has_key(primer):
                import tkMessageBox, copy
                if tkMessageBox.askyesno("Overwrite primer?",
                                         "The file contains a primer named '%s', but your primer database already\ncontains a primer with this name.\nDo you want to overwrite your existing primer with the new one?" %primer):
                    self.parent.data['primer_dict'][primer]=copy.deepcopy(pdict[primer])
            else:
                self.parent.data['primer_dict'][primer]=copy.deepcopy(pdict[primer])
        self.show_pDB_contents()
        return

    #
    # Write currently selected primers to a space-delimited file
    #
    def write_current_primers(self,filename='current_order.txt'):
        """Write currently selected primers to a space-delimited file
           This is formatted to be uploaded to the MGW website for ordering"""
        #self.save_to_text(filename='current_order.csv',separator=' ',nodesc=1)

        if not filename:
            print 'no filename given'
            return

        if not self.primers.curselection():
            print 'primer database is empty or no selection'
            return

        # Get currently selected primers into a list
        tmp_selection=self.primers.curselection()
        selection=[]
        for num in tmp_selection:
            selection.append(int(num))
        primers={}
        first=1
        # Now figure out the selected primers and store them in primers
        for primer_num in selection:
            name_selected=self.primer_order[primer_num]
            primers[name_selected]=self.parent.data['primer_dict'][name_selected]


        #fd=open(filename,'wb')

        self.DumpPDB = []
        import mutation
        for entry in primers:
            tmp = []
            tmp.append(entry+'_for')
            tmp.append(' ')
            tmp.append(primers[entry]['sequence'])
            reverse=mutation.get_reverse_complementary(primers[entry]['sequence'])
            tmp.append('\n')
            tmp.append(entry+'_rev')
            tmp.append(' ')
            tmp.append(reverse)
            tmp.append('\n')
            self.DumpPDB.append(tmp)
        try:
            fd=open(filename,'a')
            for line in self.DumpPDB:
                str = ''
                for S in line:
                    str = str+S
                fd.write(str+'\n')
            fd.close()
            print "File written ",filename
        except:
            print "error: could not write file ",filename

        from PEATDB.textFrame import textFrame
        self.order_frame = textFrame(self, 'Current primers')
        self.order_frame.load_text(self.DumpPDB)
        fd.close()

        return

    #
    # Creates current primers file and uploads to website - unused
    #
    def upload_primers(self):
        """Creates current primers file and uploads to website"""
        self.write_current_primers()
        self.get_login_info()
        #self.login_MWG()
        return

    #
    # Dialog to get login info for MWG site - unused
    #
    def get_login_info(self):
        """Dialog to get login info for MWG site"""
        self.email_var=''
        self.passwd_var=''
        self.logininfo_win=Toplevel()
        self.logininfo_win.geometry('+200+350')
        self.logininfo_win.title('Login to MWG')

        row=1
        email_lbl=Label(self.logininfo_win,text='Email:')
        email_lbl.grid(row=row,column=0,padx=3,pady=2)
        emailentry=Entry(self.logininfo_win,textvariable=self.email_var,width=50)
        emailentry.grid(row=row,column=1, sticky='news', padx=3,pady=2)
        row=2
        passwd_lbl=Label(self.logininfo_win,text='Password:')
        passwd_lbl.grid(row=row,column=0,padx=3,pady=2)
        passwdentry=Entry(self.logininfo_win,textvariable=self.passwd_var,width=50,show='*')
        passwdentry.grid(row=row,column=1, sticky='news', padx=3,pady=2)
        row=3
        close=Button(self.logininfo_win,text='Close',command=self.close_login)
        close.grid(row=row,column=0,sticky='news',padx=2,pady=2)
        ok=Button(self.logininfo_win,text='OK',command=self.login_MWG)
        ok.grid(row=row,column=1,sticky='news',padx=2,pady=2)
        return

    #
    # -------
    #
    def close_login(self,event=None):
        """Close the login info window"""
        self.logininfo_win.destroy()
        return

    #
    # Send a login attempt to the MWG server and return the html file - unused
    #
    def login_MWG(self,URL=None,event=None):
        """Send a login attempt to the MWG server and return the html file"""
        self.HOST='enzyme.ucd.ie'
        self.PORT=80
        self.URL = 'http://www.someserver.com/somepath/someprotectedpage.html'
        import urllib2
        req = urllib2.Request(self.URL)
        try:
            handle = urllib2.urlopen(req)
        except IOError, e:
            if hasattr(e, 'code'):
                if e.code != 401:
                    print 'We got a 401 error'
                    print e.code
                else:
                    print e.headers
                    print e.headers['www-authenticate']

        return




