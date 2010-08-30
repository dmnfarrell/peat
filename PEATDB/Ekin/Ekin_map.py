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
import Pmw
import os
import numpy

class Ekin_map_annotate:

    def map_datatab2structure(self):
        """If the PEATDB record has a structure, then we allow the user to map each datatab
        to a specific part of the protein.

        One can map a datatab to an atom, a residue, a chain, or define a structural group and map to it"""
        if not self.parent:
            import tkMessageBox
            tkMessageBox.showinfo("No PEAT",
                                  "This option is only available when Ekin is started from PEAT",
                                  parent=self.ekin_win)
            return
        #
        # Do we have a record name
        #
        if not self.protein:
            import tkMessageBox
            tkMessageBox.showinfo("No PEAT record",
                                  "This option is only available when Ekin has been started by clicking a PEAT record",
                                  parent=self.ekin_win)
            return
        #
        # Is there a structure?
        #
        error=None
        if not self.parent.data.has_key('DBinstance'):
            error=1
        else:
            DB=self.parent.data['DBinstance'].DB
            if not DB[self.protein].has_key('Structure'):
                error=1
            else:
                print 'Trying to get PDB'
                self.pdblines,X=self.parent.get_structure(self.protein,'Structure')
                if not self.pdblines:
                    error=1
        if error:
            import tkMessageBox
            tkMessageBox.showinfo("No Structure in PEAT",
                                  "This option is only available when the PEAT record has a structure",
                                  parent=self.ekin_win)
            return
        #
        # Open the mapping window
        #
        mapper_win=Toplevel()
        mapper_win.title('Map datatab to structure. %s - %s' %(self.protein,self.field))
        self.set_geometry(self.ekin_win,mapper_win)

        #
        # Mapping Manager
        #
        row=0
        Label(mapper_win,text='Mapping Manager',bg='lightblue').grid(row=row,column=0,columnspan=3,sticky='news')
        row=row+1
        Label(mapper_win,textvariable=self.currentdataset.get()).grid(row=row,column=0,columnspan=3,sticky='news')
        #
        # Headers
        #
        #row=row+1
        #Label(mapper_win,text='Structural group type').grid(row=row,column=0,sticky='news')
        #Label(mapper_win,text='Structural element').grid(row=row,column=1,sticky='news')
        #Label(mapper_win,text='Datatab property').grid(row=row,column=2,sticky='news')
        #
        # Structural groupings for this protein
        #
        #if not DB[self.protein].has_key('structgroups'):
        #    DB[self.protein]['structgroups']={}
        #structgroups=DB[self.protein]['structgroups'].keys()
        #
        # Load the residue definitions
        #
        import Protool.mutate
        self.M_instance=Protool.mutate.Mutate(onlydefs=1)
        self.AAdefs=self.M_instance.aadefs
        #
        # Struct group types
        #
        row=row+1
        listbox_height=5
        self.group_type_box = Pmw.ScrolledListBox(mapper_win,
                                                  items=['Residues','Atoms','Titratable groups'],
                                                  labelpos='nw',
                                                  label_text='Group type',
                                                  listbox_height = listbox_height,
                                                  usehullsize = 1,
                                                  hull_width = 200,
                                                  hull_height = 100,
                                                  selectioncommand=self.update_elements)
        self.group_type_box.grid(row=row,column=0,columnspan=1,sticky='news')
        self.group_type_box.configure(listbox_bg='white')
        self.group_type_box.configure(listbox_selectmode='single')
        self.group_type_box.configure(listbox_exportselection=0)
        #
        #
        # Dropdown list of elements of each structgroup type
        #
        self.group_elements_box = Pmw.ScrolledListBox(mapper_win,
                                              items=[],
                                              labelpos='nw',
                                              label_text='Group Elements',
                                              listbox_height = listbox_height,
                                              usehullsize = 1,
                                              hull_width = 200,
                                              hull_height = 100)
        self.group_elements_box.grid(row=row,column=1,columnspan=1,sticky='news')
        self.group_elements_box.configure(listbox_bg='white')
        self.group_elements_box.configure(listbox_selectmode='extended')
        self.group_elements_box.configure(listbox_exportselection=0)

        # Parameters that we can map to structgroups
        import Fitter
        self.FIT=Fitter.FITTER('1 pKa 2 Chemical shifts',self)

        self.dataprops=['Data source']+self.FIT.parameter_names
        self.data_prop_box = Pmw.ScrolledListBox(mapper_win,
                                              items=self.dataprops,
                                              labelpos='nw',
                                              label_text='Data properties',
                                              listbox_height = listbox_height,
                                              usehullsize = 1,
                                              hull_width = 200,
                                              hull_height = 100)
        self.data_prop_box.grid(row=row,column=2,columnspan=1,sticky='news')
        self.data_prop_box.configure(listbox_bg='white')
        self.data_prop_box.configure(listbox_selectmode='extended')
        self.data_prop_box.configure(listbox_exportselection=0)
        #
        # List of existing mappings
        #
        row=row+1
        datatab=self.currentdataset.get()
        print 'Loading this datatab in mapper',datatab
        mappings=self.get_structmappings(datatab)
        self.mapping_box = Pmw.ScrolledListBox(mapper_win,
                                               items=mappings,
                                               labelpos='nw',
                                               label_text='Existing mappings',
                                               listbox_height = 6,
                                               usehullsize = 1,
                                               hull_width = 200,
                                               hull_height = 200)
        self.mapping_box.grid(row=row,column=0,columnspan=3,sticky='news')
        self.mapping_box.configure(listbox_selectmode='single')
        self.mapping_box.configure(listbox_bg='white')
        #
        # Buttons
        #
        row=row+1
        Button(mapper_win,text='Create mapping',bg='lightgreen',borderwidth=2, relief=GROOVE, command=self.create_mapping).grid(row=row,column=0,sticky='news',padx=2,pady=2)
        Button(mapper_win,text='Delete mapping',bg='yellow',borderwidth=2, relief=GROOVE, command=self.delete_mapping).grid(row=row,column=1,sticky='news',padx=2,pady=2)
        Button(mapper_win,text='Export',bg='#CFECEC',borderwidth=2, relief=GROOVE, command=self.export_dialog).grid(row=row,column=2,sticky='news',padx=2,pady=2)

        row=row+1
        Button(mapper_win,text='Close',borderwidth=2, relief=GROOVE,command=self.close_mapper_window).grid(row=row,column=1,columnspan=2,sticky='news',padx=2,pady=2)
        #
        # Structural group manager
        #
        #row=row+1
        #Label(mapper_win,text='Structural Group Manager',bg='lightblue').grid(row=row,column=0,columnspan=3,sticky='news')
        #import os, sys
        #PEAT_dir=os.path.split(__file__)[0]
        #sys.path.append(PEAT_dir)
        #import protein_selector
        #row=row+1
        #SEL=protein_selector.select_residue(mapper_win,self.pdblines)
        #SEL.box.grid(row=row,column=0)
        ##
        #row=row+1
        #Label(mapper_win,text='Atoms').grid(row=row,column=1)
        #row=row+1
        #Button(mapper_win,text='Create new structural grouping',command=self.create_new_structgroup).grid(row=row,column=0)
        #Button(mapper_win,text='Add to structural grouping',command=self.add_to_structgroup).grid(row=row,column=1)
        #Button(mapper_win,text='Close',command=mapper_win.destroy).grid(row=row,column=2,sticky='news')
        mapper_win.rowconfigure(2,weight=1)
        self.mapper_win=mapper_win
        self.mapper_win.transient(master=self.ekin_win)
        return

    #
    # ----
    #

    def close_mapper_window(self):
        """Close the mapping window and delete references to it"""
        self.mapper_win.destroy()
        if hasattr(self,"mapper_win"):
            delattr(self,"mapper_win")
        return


    #
    # ----
    #

    def update_elements(self):
        """Insert a new dropdown list for the element"""
        #
        # Get the group type
        #
        elements=None
        group_type=self.group_type_box.getcurselection()[0]
        import Protool
        if group_type=='Residues':
            P=Protool.structureIO()
            P.parsepdb(self.pdblines)
            residues=P.residues.keys()
            residues.sort()
            elements=[]
            for res in residues:
                elements.append('%s %s' %(res,P.resname(res)))
        elif group_type=='Atoms':
            P=Protool.structureIO()
            P.parsepdb(self.pdblines)
            atoms=P.atoms.keys()
            for res in P.residues.keys():
                resname=P.resname(res)
                if self.AAdefs.has_key(resname):
                    defatoms=self.AAdefs[resname]['atoms']
                    #print defatoms
                    for defatom,coord,dummy in defatoms:
                        atom_name='%s:%s' %(res,defatom)
                        if not P.atoms.has_key(atom_name):
                            atoms.append(atom_name)
                            #print 'Adding',atom_name
            atoms.sort()
            elements=[]
            for at in atoms:
                elements.append(at)
        elif group_type=='Titratable groups':
            P=Protool.structureIO()
            P.parsepdb(self.pdblines)
            P.get_titratable_groups()
            titgrps=P.titratable_groups.keys()
            titgrps.sort()
            elements=[]
            for res in titgrps:
                for titgrp in P.titratable_groups[res]:
                    name='%s %s' %(res,titgrp['name'])
                    elements.append(name)
        else:
            print 'Unkown group type',group_type
        #
        # Make the new dropdown list
        #
        if elements:
            self.group_elements_box.setlist(elements)
        return

    #
    # -----
    #

    def create_mapping(self):
        """Create the mapping"""
        g_type=self.group_type_box.getcurselection()
        if len(g_type)==0:
            return
        g_type=g_type[0]
        g_elements=self.group_elements_box.getcurselection()
        props=self.data_prop_box.getcurselection()
        #
        if not getattr(self,'structmappings',None):
            self.structmappings={}
        datatab=self.currentdataset.get()
        if not self.structmappings.has_key(datatab):
            self.structmappings[datatab]={}
        #
        # Get the dict of current mappings
        #
        curmappings=self.structmappings[datatab]
        map_keys=curmappings.keys()
        map_keys.sort()
        #
        # Get the number of the last mapping
        #
        last_num=0
        if len(map_keys)>0:
            last_num=map_keys[-1]
        #
        # Add the new mapping
        #
        if props and g_elements and g_type:
            self.structmappings[datatab][last_num+1]={'Group type':g_type,'Group elements':g_elements,'Data property':props}
        #
        # Display the updated list of mappings
        #
        mappings=self.get_structmappings(datatab)
        self.mapping_box.setlist(mappings)
        return

    #
    # ----
    #

    def get_structmappings(self,datatab):
        """Get a printable list of structural mappings for this datatab"""
        if not getattr(self,'structmappings',None):
            return []
        if self.structmappings.has_key(datatab):
            map_keys=self.structmappings[datatab].keys()
            map_keys.sort()
            mappings=[]
            for map_key in map_keys:
                thismap=self.structmappings[datatab][map_key]
                mappings.append('%2d: %s mapped to type "%s" elements %s' %(map_key,thismap['Data property'],thismap['Group type'],thismap['Group elements']))
        else:
            mappings=[]
        return mappings

    #
    # -----
    #

    def delete_mapping(self):
        """Delete a structmapping"""
        delete=self.mapping_box.getcurselection()
        if len(delete)==0:
            print 'length is zero'
            return
        delete=str(delete[0])
        number=int(delete.split(':')[0])
        print 'NUMBER',number
        datatab=self.currentdataset.get()
        print self.structmappings.keys()
        if self.structmappings.has_key(datatab):
            if self.structmappings[datatab].has_key(number):
                del self.structmappings[datatab][number]
                mappings=self.get_structmappings(datatab)
                self.mapping_box.setlist(mappings)
        return

    #
    # -----
    #

    def update_mapping_window(self):
        """Update the mapping window when we change datatabs"""
        #
        # Update list of current mappings
        #
        datatab=self.currentdataset.get()
        mappings=self.get_structmappings(datatab)
        self.mapping_box.setlist(mappings)
        #
        # Update List of parameters
        #
        dataprops=['Data source']+self.FIT.parameter_names
        self.data_prop_box.setlist(dataprops)
        return

    def get_assigned(self):
        """Get all unique assigned elements from the mapping dict"""
        if not getattr(self,'structmappings',None):
            return []
        assigned=[]
        for key in self.structmappings.keys():
            for val in self.structmappings[key].keys():
                elements=self.structmappings[key][val]['Group elements']
                for e in elements:
                    if not e in assigned:
                        assigned.append(e)

        return assigned
    #
    # -----
    #
    def export_dialog(self):
        if hasattr(self, 'export_win'):
            if self.export_win != None :
                self.export_win.deiconify()
                return
        self.export_win=Toplevel()
        self.export_win.title('Export mappings')
        self.set_geometry(self.ekin_win,self.export_win)
        #self.setgeometry(self.ekin_win,self.export_win)
        self.grouptype = StringVar()    #group type
        grptypes=['Residues','Atoms','Titratable groups','Any']
        self.grouptype.set(grptypes[0])
        self.assignedto = StringVar()   #titratable group assigned
        #self.expdataprops=['Data source']+self.FIT.parameter_names
        self.expdataprops=['Data source','pK','span','offset']
        self.dataprop = StringVar()     #required property
        self.dataprop.set(self.expdataprops[0])
        elements=self.get_assigned()
        elements.append('All')
        elements.sort()
        self.assignedto.set(elements[0])

        row=0
        help=Label(self.export_win,text='Use the list of currently assigned mappings to select\n'
                                         +'an assigned residue/element from.\n'
                                         +'A file will be created for the chosen group element',
                                         bg='#CFECEC' )
        help.grid(row=row,column=0,columnspan=2,sticky='news',padx=2,pady=2)
        row=1
        #drop down labels for grp element, data property and assignedto
        Label(self.export_win,text='Assigned:').grid(row=row,column=0,sticky='news',padx=2,pady=2)
        w = OptionMenu(self.export_win, self.assignedto, *elements)
        w.grid(row=row,column=1,sticky='news',padx=2,pady=2)

        '''row=row+1
        Label(self.export_win,text='group type:').grid(row=row,column=0,sticky='news',padx=2,pady=2)
        w = OptionMenu(self.export_win, self.grouptype, *grptypes)
        w.grid(row=row,column=1,sticky='news',padx=2,pady=2)'''

        row=row+1
        Label(self.export_win,text='data property:').grid(row=row,column=0,sticky='news',padx=2,pady=2)
        print self.dataprops
        w = OptionMenu(self.export_win, self.dataprop, *self.expdataprops)
        w.grid(row=row,column=1,sticky='news',padx=2,pady=2)

        row=row+1
        Button(self.export_win,text='Cancel',bg='#CFECEC',borderwidth=2, relief=GROOVE, width=10,
                command=self.close_exp_dialog).grid(row=row,column=0,sticky='news',padx=2,pady=2)
        Button(self.export_win,text='Go',bg='#CFECEC',borderwidth=2, relief=GROOVE, width=10,
                command=self.export_as_csv).grid(row=row,column=1,sticky='news',padx=2,pady=2)

        return

    def close_exp_dialog(self):
        if hasattr(self,'export_win'):
            self.export_win.destroy()
            self.export_win=None
        return

    def choose_savedir(self):
        """Get a directory to save to"""
        import tkFileDialog, os
        if self.defaultsavedir == None:
            self.defaultsavedir = os.getcwd()
        dirname=tkFileDialog.askdirectory(parent=self.export_win,
                                          initialdir=self.defaultsavedir)
        if not dirname:
            print 'Returning'
            return NoneType

        return dirname
    #
    # -----
    #
    def export_as_csv(self):
        """export struct mapping for specific filters as csv"""
        #prompt user for save dir
        savedir = self.choose_savedir()
        if savedir==None:
            return
        if self.currplatform == 'Windows':
            print 'using windows'

        import List_Utils
        #sub function for tidiness
        def getexplist(assignedto):

            reslist={}
            reskeys=[]
            for key in self.structmappings.keys():
                for n in self.structmappings[key].keys():
                    #check if any dataprop list element contains the key eg 'pK' in pK1, pK2 etc..
                    datapropkey = List_Utils.elements_contain(self.structmappings[key][n]['Data property'], self.dataprop.get())
                    if datapropkey != None:
                        #try to extract the value from the ekin dataset
                        val = self.get_dataprop_value(key, datapropkey)
                        print 'found ',val,' for ', datapropkey
                        #print 'val: ', val
                        #iterate over group elements list
                        elements=self.structmappings[key][n]['Group elements']
                        for e in elements:
                            if assignedto in e:
                                reslist[key] = ([key,val])
                                reskeys.append(key)

            if len(reslist.keys())==0:
                return
            #write the list to a csv file, first add heading
            import string
            #remove whitespace
            name=string.join(assignedto.split(), '')
            name=name.replace(':', '')
            if self.currplatform == 'Windows':
                filename = savedir+'/'+name+'.csv'
            else:
                filename = os.path.join(savedir, name+'.csv')
            print filename
            writer = open(filename, "wb")
            writer.write(assignedto+'\n')
            import csv
            csvwriter = csv.writer(open(filename, "a"))
            keyssorted = self.sort_by_Num(reskeys)
            #print reslist
            #print keyssorted
            p=[];names=[]
            #use key sorted mapping to list residues by number
            for item in keyssorted:
                k=item[1]
                csvwriter.writerow(reslist[k])
                p.append(reslist[k][1])
                names.append(k)
            writer.close()
            #do a plot and save to same dir as file
            try:
                import pylab
            except:
                return
            f=pylab.figure(figsize=(10,4))
            pylab.rc("font", family='serif')
            a=f.add_subplot(111)
            ind=numpy.arange(len(names))
            a.bar(ind, p , linewidth=0.5)
            a.set_xticks(ind)
            a.set_ylabel(self.dataprop.get())
            a.set_title(name+' assignments')
            a.set_xticklabels(names, rotation='vertical', size=5)
            f.savefig(savedir+'/'+name+'.png',dpi=300)
            return

        if self.assignedto.get() == 'All':
            for a in self.get_assigned():
                getexplist(a)
        else:
            getexplist(self.assignedto.get())

        self.close_exp_dialog()
        return

    #
    # -----
    #
    def get_dataprop_value(self, key, dataprop):
        """Annoying but necessary helper func to get value of assigned property
           from the ekin fit data"""
        tabnum = self.currentdataset.get()
        if self.fitter_data.has_key(key):
            fitdata = self.fitter_data[key]
        else:
            return None
        model = fitdata['model']
        #extracts index number from fit model field name
        i = self.FIT.get_param_index(dataprop, model)
        print tabnum, key
        print fitdata, dataprop, i
        if i!=None:
            val = fitdata[i]
            return val

    #
    # -----
    #
    def create_new_structgroup(self):
        return

    #
    # ------
    #

    def add_to_structgroup(self):
        return



    def sort_by_Num(self, p):
        """Sort text keys by contained numbers - should be put in utils class"""
        splitkeys={}
        import re
        r=re.compile('\D')
        for k in p:
            splitkeys[k]=int(r.split(k)[1])
        items = splitkeys.items()
        items = [(v, k) for (k, v) in items]
        items.sort()
        return items

