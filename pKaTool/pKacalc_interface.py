#!/usr/bin/env python
#
# pKaTool - analysis of systems of titratable groups
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

from Tkinter import *
import Pmw
import tkMessageBox, tkFileDialog

class Main:

    def load_pdb(self):
        """Open a window to load a PDB file with the aim of calculating pKa values

        Get the PDB file, and see if a calculation has been carried out"""
        import os

        if not self.pdbfile:
            pdbfilename=tkFileDialog.askopenfilename(defaultextension='.pdb',
                                                     initialdir=os.getcwd(),
                                                     filetypes=[("PDB file","*.pdb"),
                                                                ("All files","*.*")],
                                                     parent=self.master)
            if not pdbfilename:
                return
        else:
            pdbfilename=self.pdbfile
        #
        #
        #
        self.pdbfilename=os.path.join(os.getcwd(),pdbfilename)
        if not os.path.isfile(self.pdbfilename):
            tkMessageBox.showwarning("File not found",
                                     'File not found: %s' %pdbfilename
                                     ,parent=self.master)
            return
        #
        # Set the name for it
        #
        self.label_count=0
        self.labels='ABCDEFGHIJKLMNOP'
        label=self.labels[self.label_count]
        #
        # Load the PDB file in Protool
        #
        import Protool
        self.Protool_instance=Protool.structureIO()
        self.Protool_instance.readpdb(self.pdbfilename)
        #
        # Setup the calculation arrays
        #
        import pKaIO
        self.calcs[self.pdbfilename]={'instance':pKaIO.pKaIO(self.pdbfilename),'selected_groups':[],'label':label}
        self.label_count=self.label_count+1
        self.X=self.calcs[self.pdbfilename]['instance']
        #
        # Open the file and assess status
        #
        if self.X.assess_status():
            import tkMessageBox
            answer=tkMessageBox.askyesno('WHAT IF pKa calculation found',
                                          'A WHAT IF pKa calculation has been found for this PDB file.\nDo you want to load the results from this calculation?',
                                          parent=self.master)
            if answer:
                self.load_titration()
            else:
                self.pdbfilename_var.set('PDB file: %s' %os.path.split(pdbfilename)[1])
                self.pka_status_var.set('pKa calc status: No calculation performed')
                del self.calcs[pdbfilename]
        else:
            self.pdbfilename_var.set('PDB file: %s' %os.path.split(pdbfilename)[1])
            self.pka_status_var.set('pKa calc status: No calculation performed')
            del self.calcs[pdbfilename]
        return

    #
    # ----
    #

    def full_pKa_updater(self,message):
        self.status2.set(message)
        self.master.update_idletasks()
        return

    def pKaTool_pKa(self):
        """Perform a dirt-cheap pKa calculation using Jens' simple, simple method"""
        self.pHstart=0.1
        self.pHstop=20.0
        self.pHstep=0.1
        self.chwin.destroy()
        self.simple_win=Toplevel()
        self.simple_win.transient(self.master)
        self.simple_win.title('pKaTool cheap pKa calculation')
        self.set_geometry(self.master,self.simple_win)
        self.status=StringVar()
        self.status.set('Status: Initialising pKa calculation routines')
        Label(self.simple_win,textvariable=self.status).grid(row=5,column=0,columnspan=2)
        self.status2=StringVar()
        self.status2.set('---')
        Label(self.simple_win,textvariable=self.status2,bg='white').grid(row=6,column=0,columnspan=2)
        self.update()
        #
        # Do we have a PDB file?
        #
        if not self.pdbfilename:
            tkMessageBox.showwarning("No PDB file",
            "You did not load a PDB file",
            parent=self.simple_win)
            self.simple_win.destroy()
            return
        #
        # Start the pKa calculation
        #
        import Protool
        import Protool.errors
        P=Protool.structureIO()
        self.status.set('Reading PDB file')
        self.update()
        P.readpdb(self.pdbfilename)
        self.status.set('Removing non AA atoms')
        P.Remove_All_NonAminoAcids()
        try:
            self.status.set('Calculating site-site interaction energies')
            self.update()
            matrix=P.calculate_matrix(self.pdbfilename,updater=self.full_pKa_updater)
            self.status.set('Calculating desolvation energies')
            self.update()
            desolv=P.calculate_desolvation(self.pdbfilename,updater=self.full_pKa_updater)
            #
            self.status.set('Calculating background interaction energies')
            self.update()
            backgr=P.calculate_background(self.pdbfilename,updater=self.full_pKa_updater)
        except Protool.errors.AtomNotFoundError,inst:
            tkMessageBox.showwarning("Atom missing",
                                     "At least one atom is missing in the PDB file.\nMake sure that all atoms are present before starting a pKa calculation.\nYou can e.g. use the pdb2pqr server for that.\nError:\n%s" %inst,
                                     parent=self.simple_win)
            self.simple_win.destroy()
            return

        self.status.set('Calculating titration using Monte Carlo sampling')
        self.update()

        import pKa_calc
        try:
            C=pKa_calc.Monte_Carlo_CPP()
            C.test()
        except:
            tkMessageBox.showwarning("C++ module not found",
                                     'C++ module pMC not found. I will revert to a python implementation of the Tanford-Roxby algorithm.\n\nIf you are running linux, then you should try to cd to the pKaTool directory and type "make"',
                                     parent=self.master)
            self.status.set('Calculating titration curves using Tanford-Roxby algorithm')
            self.update()
            C=pKa_calc.Tanford_Roxby()
        #
        # Load the values
        #
        C.matrix=matrix
        C.desolv=desolv
        C.backgr=backgr
        C.calc_pKas(mcsteps=200000,phstep=self.pHstep,phstart=self.pHstart,phend=self.pHstop,complete_pka=1)
        #
        # Write the titration curve file and the PKA.DAT file
        #
        import pKaIO
        X=pKaIO.pKaIO()
        X.write_titration_curve(self.pdbfilename+'.TITCURV.DAT',C.prot_states)
        X.write_pka(self.pdbfilename+'.PKA.DAT',C.pka)

        self.status.set('Done')
        self.update()
        #
        # Load the results
        #
        self.simple_win.destroy()
        self.chwin.destroy()
        self.load_titration(pdbfilename=self.pdbfilename)
        tkMessageBox.showwarning("pKa calculation done",'The pKa calculation has finished',parent=self.master)
        return

    #
    # ----
    #

    def WI_pKa(self):
        """Setup for a WHAT IF pKa calculation"""
        tkMessageBox.showwarning("Not implemented",'Not implemented yet',parent=self.calc_window)
        return

    #
    # -----
    #

    def load_frompKD(self,event=None):
        """Load calculated pKa values from the pKD server"""
        self.pKD_IFwin=Toplevel()
        self.set_position(self.pKD_IFwin)
        self.pKD_IFwin.title('pKD server connection')
        self.URL='/cgi-bin/pKD/pKaDesign_server.py'
        self.HOST='enzyme.ucd.ie'
        self.PORT='80'
        Label(self.pKD_IFwin,text='pKD server at: http://%s:%s%s' %(self.HOST,self.PORT,self.URL),bg='green').grid(row=0,column=0)
        #
        # Get the list of proteins that are prepared
        #
        request ={}
        X=HTTP_handler(self.URL,self.HOST,self.PORT)
        request['action']='get_prepared_proteins'
        request['client']='PKATOOL'
        X.load(request)
        lines=X.request(return_type='file-object')
        #
        # Parse the XML
        #
        proteins=[]
        from xml.dom import minidom
        xmldoc=minidom.parse(lines)
        prot_sec=xmldoc.firstChild
        for calc in prot_sec.childNodes:
            ID=calc.attributes['id'].value
            for attr in calc.childNodes:
                if attr.nodeName=='info':
                    import string
                    if len(attr.childNodes)>0:
                        info=string.strip(attr.childNodes[0].data)
                    else:
                        info='No info'
            proteins.append([ID,info,None])
        import string
        #
        # Calculations found
        #
        Label(self.pKD_IFwin,text='Calculations on the pKD server').grid(row=1,column=0)
        row=2
        yscrollbar=Scrollbar(self.pKD_IFwin,orient='vertical',width=14)
        yscrollbar.grid(row=row,column=2,sticky='nws',rowspan=5)
        self.calculations=Listbox(self.pKD_IFwin,
                                  bg='white',
                                  fg='black',
                                  height=10,width=15,yscrollcommand=yscrollbar.set,
                                  selectmode=MULTIPLE)
        self.calculations.grid(row=row,column=0,columnspan=2,sticky='news',padx=2,rowspan=5)
        yscrollbar.config(command=self.calculations.yview)
        #
        # Sort the calcs and insert them
        #
        self.calculations_ready=proteins
        self.calculations_ready.sort()
        self.calcs_displayed=self.calculations_ready
        for calcname,info,setup in self.calculations_ready:
            self.calculations.insert(END,'%6s: %25s' %(calcname,info[:45]))
        #
        # Box with details
        #
        self.searchtext=StringVar()
        l=Label(self.pKD_IFwin,text='Search ')
        l.grid(row=row,column=3,sticky='nes')
        self.searchbox=Entry(self.pKD_IFwin,textvariable=self.searchtext)
        self.searchbox.grid(row=row,column=4,sticky='w')
        self.searchbox.bind('<Return>',self.do_search)
        #
        #
        row=row+1
        lbl=Label(self.pKD_IFwin,text='this is where details on the calcs will go in the future')
        lbl.grid(row=row,column=3,columnspan=5)
        #
        # Buttons for selecting
        #
        row=row+6
        Button(self.pKD_IFwin,text='Select All',command=self.pKDselect_all).grid(row=row,column=0)
        Button(self.pKD_IFwin,text='Clear selection',command=self.pKDselect_none).grid(row=row,column=1)
        #
        # Button for loading calculation
        #
        row=row+1
        self.singlebutton=Button(self.pKD_IFwin,text='Load calculation(s)',command=self.do_load_pKDcalcs)
        self.singlebutton.grid(row=row,column=0)
        Button(self.pKD_IFwin,text='Cancel',command=self.pKD_IFwin.destroy).grid(row=row,column=1)
        return

    #
    # ----
    #

    def pKDselect_all(self):
        """Select all of the calcs displayed in the listbox"""
        self.calculations.selection_set(0,END)
        return

    #
    # -----
    #

    def pKDselect_none(self):
        """Clear the selection"""
        self.calculations.selection_clear(0,END)
        return

    #
    # ----
    #

    def do_search(self,event=None):
        """Search for the calcs where a certain text is found"""
        self.calculations.delete(0,END)
        self.calcs_displayed=[]
        text=self.searchtext.get()
        for calcname,info,setup in self.calculations_ready:
            if calcname.find(text)!=-1 or info.find(text)!=-1:
                self.calcs_displayed.append([calcname,info,setup])
                self.calculations.insert(END,'%6s: %25s' %(calcname,info[:45]))
        return

    #
    # ----
    #

    def get_pKD_response(self,command,URL=None):
        """Send a command to the pKD server and return the response"""
        request ={}
        for c in command.keys():
            request[c]=command[c]
        request['client']='PKATOOL'
        #
        if not URL:
            URL=self.URL
        X=HTTP_handler(URL,self.HOST,self.PORT)
        X.load(request)
        #
        lines=X.request().split('\n')
        response=[]
        import string
        for count in range(len(lines)):
            line=lines[count]
            if string.strip(line)=='PKATOOLOUTPUT:':
                response=lines[count+1:]
        return response

    #
    # ----
    #


    #
    # ----
    #

    def do_load_pKDcalcs(self):
        """Load pKa calculations from the pKD server into pKaTool"""
        selection=self.calculations.curselection()
        if len(selection)==0:
            return
        self.pKDname=StringVar()
        self.pKDname.set('Ready')
        Label(self.pKD_IFwin,textvariable=self.pKDname).grid(row=26,column=0)
        count=0
        for sel in selection:
            count=count+1
            self.do_load_singlepKDcalc(int(sel),count,len(selection))
        self.pKD_IFwin.destroy()
        return

    #
    # ----
    #

    def do_load_singlepKDcalc(self,selection,count=1,total=1):
        """Load a single pKD calculation"""
        print 'Loading zipfile'
        calc_selected,info,setup=self.calcs_displayed[selection]
        self.pKDname.set('Loading: %s.... (%3d of %3d)' %(calc_selected,count,total))
        self.pKD_IFwin.update_idletasks()
        #
        # Get the zipfile from the server
        #
        zip_URL=self.get_pKD_response({'action':'get_pka_calc','PDBID':calc_selected})[0]
        #
        import urllib, zipfile, tempfile, os
        f = urllib.urlopen(zip_URL)
        tmp=tempfile.TemporaryFile()
        tmp.writelines(f.read())
        f.close()
        zfile=zipfile.ZipFile(tmp,'r')
        tmpdir=tempfile.mkdtemp()
        pdbfilename=None
        for file in zfile.namelist():
            newname=os.path.join(tmpdir,file)
            fd=open(newname,'w')
            fd.writelines(zfile.read(file))
            fd.close()
            if newname[-4:]=='.pdb':
                pdbfilename=newname
        self.load_titration(pdbfilename)
        return
        
    #
    # ----
    #

    def do_load_singlepKDcalcXML(self,selection,count=1,total=1):
        print 'Loading XMLfile'
        calc_selected,info,setup=self.calcs_displayed[selection]
        self.pKDname.set('Loading: %s.... (%3d of %3d)' %(calc_selected,count,total))
        self.pKD_IFwin.update_idletasks()
        #
        # Get the XML data from the server
        #
        request ={}
        X=HTTP_handler(self.URL,self.HOST,self.PORT)
        request['action']='get_xml_pKa_calc'
        request['client']='PKATOOL'
        request['PDBID']=calc_selected
        X.load(request)
        lines=X.request(return_type='file-object')
        #print lines.read()
        #
        # Parse the XML
        #
        proteins=[]
        from xml.dom import minidom
        xmldoc=minidom.parse(lines)
        pKas=xmldoc.firstChild
        for calc in pKas.childNodes:
            ID=calc.attributes['id'].value
            for attr in calc.childNodes:
                pKaval=attr
            print '----'
    
        return


    #
    # -----
    #

    def get_titratable_groups(self):
        """Get the titratable groups in the current file, with the selected method"""
        engine=self.calc_engine.get()
        if engine=='PDB2pKa':
            self.get_pdb2pKa_groups()
        elif engine=='pKaTool' or engine=='Manual energies':
            self.Protool_instance.get_titratable_groups()
            PT_grps=self.Protool_instance.titratable_groups.keys()
            PT_grps.sort()
            self.groups=[]
            for residue in PT_grps:
                for group in self.Protool_instance.titratable_groups[residue]:
                    self.groups.append(self.Protool_instance.get_titgroup_name(residue,group))
        else:
            print 'Not implemented yet'
        return

    #
    # ----
    #

    def file_manager_dialog(self):
        """Get the titratable groups as defined by the PDB2pKa routines"""
        self.file_win=Toplevel()
        self.file_win.geometry('+400+500')
        if self.pdbfilename:
            self.files_included=[[self.pdbfilename,'Protein']]
        self.files_included=[]
        self.print_included()
        return

    #
    # -----
    #

    def print_included(self):
        """Print the files that we are considering"""
        #
        # Do we have any files?
        #
        if len(self.files_included)<1:
            return
        Label(self.file_win,text='Files holding structure data').grid(row=0,column=0,columnspan=3)
        #
        # Print the file(s) that we use
        #
        import os
        row=1
        count=1
        for name,description in self.files_included:
            Label(self.file_win,text='%d' %count,bg='white').grid(row=row,column=0,sticky='news')
            Label(self.file_win,text='%s' %os.path.split(name)[1],bg='white').grid(row=row,column=1,sticky='news')
            Label(self.file_win,text='%s' %description,bg='white').grid(row=row,column=2,sticky='news')
            count=count+1
            row=row+1
        Button(self.file_win,text='Add another file',command=self.add_file).grid(row=row,column=0,columnspan=1,sticky='nws')
        Button(self.file_win,text='Done',command=self.really_getgroups).grid(row=row,column=2,columnspan=1,sticky='nes')
        return

    #
    # -----
    #

    def really_getgroups(self):
        """Divide the files listed into protein and ligand, and initialise the routines"""
        proteins=[]
        ligands=[]
        for file,ftype in self.files_included:
            if ftype=='Protein':
                proteins.append(file)
            else:
                ligands.append(file)
        #
        # Right now we can only deal with one file of each type
        #
        pdbfile=proteins[0]
        if len(ligands)>0:
            ligandfile=ligands[0]
        else:
            ligandfile=None
        #
        # Now get the groups
        #
        if self.pdb2pKa_path.get()=="":
            tkMessageBox.showerror('PDB2pKa missing',
                                   'You did not specify the path to PDB2pKa')
            return
        import sys
        sys.path.append('/home/nielsen/lib/pdb2pqr_develop/pKa')
        import pka
        myProtein, myRoutines, myForcefield,apbs_setup, ligand_titratable_groups=pka.pre_init(pdbfile=pdbfile,ff='parse',ligand=ligandfile)
        print 'Return values from pdb2pqr pka'
        print myProtein
        print myRoutines
        print myForcefield
        print apbs_setup
        print ligand_titratble_groups
        groups=pka.pKaRoutines(myProtein, myRoutines, myForcefield,apbs_setup)
        print '-------'
        print groups
        return

    #
    # ----
    #

    def add_file(self,event=None):
        """Browse for another file to add to self.included_files"""
        import os
        newfilename=tkFileDialog.askopenfilename(defaultextension='.mol2',
                                                 initialdir=os.getcwd(),
                                                 filetypes=[("Mol2 file",'*.mol2'),
                                                            ("PDB file","*.pdb"),
                                                            ("All files","*.*")],
                                                 parent=self.groupwin)
        if not newfilename:
            return

        ftype='Protein'
        if os.path.isfile(newfilename):
            if newfilename[-4:].lower()=='mol2':
                ftype='Ligand'
            #
            self.files_included.append([newfilename,ftype])
            self.print_included()
        return

    #
    # -----
    #

    def init_calc(self):
        """Initialise a pKa calculation

        This involves getting the titratable groups"""
        #
        # Do we have a PDB file?
        #
        #if not self.pdbfilename:
        #    import tkMessageBox
        #    tkMessageBox.showwarning('No PDB file',
        #                             'You have to load a PDB file before performing a pKa calculation')
        #    return
        #
        # Open the file manager dialog
        #
        self.file_manager_dialog()
        #
        # Open the window
        #
        self.initwin=Toplevel()
        self.initwin.transient(self.master)
        self.set_geometry(self.master,self.initwin)
        #
        #
        # Start window
        #
        self.initwin.title('Initialise pKa calculation')
        Label(self.initwin,text='Select the titratable groups you want to include in the calculation').grid(row=0,column=0,columnspan=5)
        Label(self.initwin,text='You cannot add or delete groups from this calculation once you click "Select these groups"').grid(row=1,column=0,columnspan=5)
        #
        # Listbox for all titratable groups
        #
        yscrollbar=Scrollbar(self.initwin,orient='vertical',width=10)
        yscrollbar.grid(row=2,column=2,rowspan=10,sticky='nws')
        height=10
        self.all_groups=Listbox(self.initwin,
                                bg='white',
                                fg='black',
                                height=height,width=15,
                                yscrollcommand=yscrollbar.set,
                                selectmode=SINGLE)
        self.all_groups.grid(row=2,column=0,columnspan=2,rowspan=10,sticky='news')
        yscrollbar.config(command=self.all_groups.yview)
        #
        # Listbox for titratable groups to be included in calculation
        #
        yscrollbar2=Scrollbar(self.initwin,orient='vertical',width=10)
        yscrollbar2.grid(row=2,column=6,rowspan=10,sticky='nws')
        height=10
        self.groups_selected=Listbox(self.initwin,
                                bg='white',
                                fg='black',
                                height=height,width=15,
                                yscrollcommand=yscrollbar2.set,
                                selectmode=SINGLE)
        self.groups_selected.grid(row=2,column=4,columnspan=2,rowspan=10,sticky='news')
        yscrollbar2.config(command=self.groups_selected.yview)
        #
        # Calculation engine setup
        #
        self.set_enecalc_engine(win=self.initwin,row=2,column=3,command=self.update_groups_selected)
        self.update_groups_selected()
        Button(self.initwin,text='Add -->',command=self.select_group).grid(row=5,column=4)
        Button(self.initwin,text='<-- Remove',command=self.remove_group).grid(row=6,column=4)
        #
        # pH start, stop and step
        #
        self.pHstart=0.1
        self.pHstop=12.0
        self.pHstep=0.1
        #
        # OK and Cancel buttons
        #
        Button(self.initwin,text='Select these groups',command=self.setup_calc).grid(row=13,column=0,columnspan=2,sticky='news')
        Button(self.initwin,text='Select all groups',command=self.setup_calc).grid(row=13,column=2,columnspan=2,sticky='news')
        Button(self.initwin,text='Add structure file(PDB/mol2)').grid(row=13,column=4)
        Button(self.initwin,text='Add extra titgroup in loaded structure')
        Button(self.initwin,text='Cancel',command=self.initwin.destroy).grid(row=13,column=5,sticky='news')
        self.status=StringVar()
        self.status.set('Waiting...')
        Label(self.initwin,text='Status:').grid(row=14,column=0)
        Label(self.initwin,textvariable=self.status).grid(row=14,column=1)
        #
        # Button for defining titratable groups
        #
        Button(self.initwin,text='Edit titgroup definitions').grid(row=14,column=0)

        return

    #
    # ----
    #

    def update_groups_selected(self):
        """Update the list of groups selected"""
        self.get_titratable_groups()
        self.all_groups.delete(0, END)
        for group in self.groups:
            self.all_groups.insert(END,group)
        return

    #
    # ----
    #

    def setup_calc(self):
        """Set up a new calculation with the present titratable groups

        Set all energies to zero
        """
        import os
        pdbfilename_short=os.path.split(self.pdbfilename)[1]
        label=self.labels[self.label_count]
        import pKaIO
        self.calcs[pdbfilename_short]={'instance':pKaIO.pKaIO(self.pdbfilename),'selected_groups':[],'label':label}
        self.label_count=self.label_count+1
        #
        # Matrix
        #
        self.status.set('Constructing matrix')
        self.master.update_idletasks()
        self.construct_empty_matrix()
        #self.calcs[pdbfilename_short]['titcurv']=self.X2.readtitcurv()
        #
        # read the pKa values
        #
        #self.calcs[pdbfilename_short]['pkavals']=self.X2.readpka()
        #
        # Store the PDB file
        #
        fd=open(self.pdbfilename)
        self.calcs[pdbfilename_short]['pdblines']=fd.readlines()
        fd.close()
        #
        # Destroy the init win
        #
        self.initwin.destroy()
        return

    #
    # ----
    #

    def construct_empty_matrix(self):
        """Construct an empty matrix and save it"""
        self.matrix={}
        for group in self.groups:
            if not self.matrix.has_key(group):
                self.matrix[group]={}
            for group2 in self.groups:
                self.matrix[group][group2]=[0.0,0.0,0.0,0.0]
        #
        # Save it
        #
        import os
        pdbfilename_short=os.path.split(self.pdbfilename)[1]
        pKa_instance=self.calcs[pdbfilename_short]['instance']
        pKa_instance.matrix=self.matrix.copy()
        pKa_instance.write_matrix(pKa_instance.matrix_file)
        return

    #
    # -----
    #
    def do_pka_calculation(self):
        """Run a full pKa calculation on the selected groups"""
        self.chwin=Toplevel()
        self.chwin.transient(self.master)
        self.set_geometry(self.master,self.chwin)
        self.chwin.title('Choose calculation type')
        self.set_enecalc_engine(self.chwin,row=0,column=0)
        Button(self.chwin,text='Run calculation',command=self.pKacalc_driver).grid(row=5,column=0)
        return

    #
    # ----
    #

    def pKacalc_driver(self):
        """Start the pKa calculation"""
        engine=self.calc_engine.get()
        if engine=='PDB2pKa':
            pass
        elif engine=='pKaTool':
            self.pKaTool_pKa()
        else:
            print 'Not implemented yet'
        return

    #
    # ----
    #

    def pdb2pKa(self):
        """Interface to pdb2pka"""
        return

    #
    # ----
    #

    def propka(self):
        """Run a pKa calculation on the propka server"""

        return

    #
    # ----
    #

    def set_enecalc_engine(self,win,row=0,column=0,command=None):
        """Set the calculation engine"""
        engines=['PDB2pKa','PropKa','WHAT IF','pKaTool','Manual energies']
        active_engines=['pKaTool','Manual energies']
        engines.sort()
        Label(win,text='Select calculation engine').grid(row=row,column=column)
        if not command:
            command=self.dummy
        for engine in engines:
            if engine in active_engines:
                Radiobutton(win,variable=self.calc_engine,value=engine,text=engine,
                            activebackground='red',command=command).grid(row=row+1,column=column,sticky='nws')
                row=row+1
        return

    #
    # ----
    #

    def dummy(self):
        """Dummy function"""
        return

    #
    # ----
    #

    def set_titcurv_method(self):
        return

    #
    # ----
    #

    def calculate_matrix_dialog(self):
        """Open a dialog for calculating or entering energies for the interaction energy matrix"""
        if not self.pdbfilename:
            tkMessageBox.showwarning('No PDB file',
                                     'You must load a PDB file before you can calculate the matrix')
            return
        #
        # Open the window
        #
        self.mwin=Toplevel()
        self.mwin.title('Interaction energy matrix')
        self.set_geometry(self.master,self.mwin)
        #
        # Canvas for the matrix
        #
        self.matrix_canvas=Pmw.ScrolledCanvas(self.mwin,
                                              borderframe = 1,
                                              labelpos = 'n',
                                              label_text = 'Interaction energy matrix',
                                              usehullsize = 1,
                                              hull_width = 800,
                                              hull_height = 800,
                                              hscrollmode='dynamic',
                                              vscrollmode='dynamic'
        )
        self.matrix_canvas.interior().configure(bg='white')
        self.matrix_canvas.grid(row=0,column=0,rowspan=10,columnspan=10,sticky='news')
        #
        # Buttons for the calculation engine
        #
        self.set_enecalc_engine(win=self.mwin,row=11,column=0)
        #
        # What do we display?
        #
        row=11
        Label(self.mwin,text='Display').grid(row=row,column=3)
        self.m_display=StringVar()
        for disp in ['Interaction energy (kT/e)','Interaction energy (dpKa)','Distance (A)','eps(effective)']:
            Radiobutton(self.mwin,variable=self.m_display,value=disp,text=disp,
                        activebackground='red',command=self.update_matrix_display).grid(row=row+1,column=3,sticky='nws')
            row=row+1
        self.m_display.set('Interaction energy (kT/e)')
        #
        # Add action buttons
        #
        row=row+1
        Button(self.mwin,text='Close',command=self.mwin.destroy).grid(row=row,column=0)
        Button(self.mwin,text='Recalculate matrix',command=self.calculate_matrix_driver).grid(row=row,column=1)
        #
        # Status field
        #
        self.status=StringVar()
        self.status.set('Waiting...')
        Label(self.mwin,text='Status:',bg='white').grid(row=14,column=0,sticky='news')
        Label(self.mwin,textvariable=self.status,bg='white').grid(row=14,column=1,sticky='news')
        #
        # Update the display
        #
        self.update_matrix_display()
        #
        # Make sure that the sliders are in place
        #
        self.matrix_canvas.resizescrollregion()
        return

    #
    # ----
    #

    def calculate_matrix_driver(self):
        """Calculate the matrix"""

        engine=self.calc_engine.get()
        if engine=='PDB2pKa':
            pass
        elif engine=='pKaTool':
            self.matrix=self.Protool_instance.calculate_matrix(filename=self.pdbfilename,updater=self.update_matrix_status)
            self.update_matrix_display()
        else:
            print 'Not implemented yet'
        return

    #
    # ----
    #

    def update_matrix_status(self,message):
        """Update the matrix calc status"""
        self.status.set(message)
        self.master.update_idletasks()
        self.mwin.update()
        return

    #
    # ----
    #

    def update_matrix_display(self):
        """Update the matrix on the screen when we choose a new calculation method"""
        #
        # Delete all object on the canvas
        #
        objects=self.matrix_canvas.find_all()
        for obj in objects:
            self.matrix_canvas.delete(obj)
        #
        # Precompute the group numbers
        #
        self.group_numbers={}
        count=0
        for group in self.groups:
            self.group_numbers[group]=count
            count=count+1
        #
        # Get the titratable groups
        #
        matrix_balloon=Pmw.Balloon(self.mwin)
        import math
        for group in self.groups:
            #
            # Labels
            #
            x,y=self.get_matrix_cell_position(group,group)
            self.matrix_canvas.create_text(x,self.matrix_y_add,text=group,anchor='s')
            self.matrix_canvas.create_text(0,y,text=group,anchor='w')
            #
            # Insert the energies
            #
            for partner in self.groups:
                x,y=self.get_matrix_cell_position(group,partner,type='box')
                if group==partner:
                    self.matrix_canvas.create_rectangle(x,y,x+self.matrix_x_add,y+self.matrix_y_add,fill='black')
                else:
                    record=self.matrix[group][partner]
                    disp=self.m_display.get()
                    intene=record[0]-record[1]-record[2]+record[3]
                    if disp=='Interaction energy (kT/e)':
                        pass
                    elif disp=='Interaction energy (dpKa)':
                        intene=intene/math.log(10)
                    elif disp=='Distance (A)':
                        import tkMessageBox
                        tkMessageBox.showinfo('Not done yet',
                                              'Be nice to Jens to get him to do this',
                                              parent=self.mwin)
                        return
                    elif disp=='eps(effective)':
                        import tkMessageBox
                        tkMessageBox.showinfo('Not done yet',
                                              'Be nice to Jens to get him to do this',
                                              parent=self.mwin)
                        return
                    if intene>0.0:
                        color='red'
                    else:
                        color='blue'
                    self.matrix_canvas.create_rectangle(x,y,x+self.matrix_x_add,y+self.matrix_y_add)
                    x,y=self.get_matrix_cell_position(group,partner)
                    handle=self.matrix_canvas.create_text(x,y,text='%5.3f' %intene,anchor='center',fill=color)
                    matrix_balloon.tagbind(self.matrix_canvas,handle,'%s-%s' %(group,partner))
                    #self.Entry=Entry(self.mwin,width=10,textvariable=energy,bg='white')
                    #self.win=self.matrix_canvas.create_window(x,y,window=self.Entry,anchor='nw')

        return

    #
    # ----
    #

    def get_matrix_cell_position(self,group,partner,type='text'):
        """Calculate the x and y coordinates the upper left corner of a matrix cell

        setting type to text gives the center of the cell
        """
        self.matrix_x_add=80
        self.matrix_y_add=25
        x=self.matrix_x_add+self.group_numbers[partner]*self.matrix_x_add
        y=self.matrix_y_add+self.group_numbers[group]*self.matrix_y_add
        if type=='text':
            y=y+int(0.5*self.matrix_y_add)
            x=x+int(0.5*self.matrix_x_add)
        return x,y



class HTTP_handler:
    """Send an HTTP request to a remote URL"""


    def __init__(self, uri, host, port):
        #
	# Initialise variables
	#
        #
        import time, os, random
        self.uri = uri
        self.host = host
        self.port = port
        self.queryString = None
        self.boundary= '%s%s_%s_%s' % \
                        ('-----', int(time.time()), os.getpid(), random.randint(1,10000))
        return


    #
    # -----
    #

    def load(self, args, headerDict=None):
        #
        # Loads multiple files
	#
        import string
        total = []
        for (name,values) in args.items():
            data = values
            hdr = []; part = []
            hdr.append('Content-Disposition: form-data; name="'+name+'"')
            part.append("%s\n\n%s" % (string.joinfields(hdr,'\n'), data))
            #print string.joinfields(hdr,'\n')
            total.append('--%s\n' % self.boundary)
            total.append(string.joinfields(part, "\n--%s\n" % self.boundary))
            total.append('\n')
        self.queryString = string.joinfields(total, '')
        return

    #
    # -----
    #

    def request(self,return_type='text'):
        """Send the request"""
        import httplib
        query = self.queryString
        contentType = 'multipart/form-data; boundary=%s' % self.boundary
        contentLength = str(len(query))
        h = httplib.HTTP()
	h.connect(self.host, self.port)
        h.putrequest('POST', self.uri)
        h.putheader('Accept', '*/*')
        h.putheader('Proxy-Connection', 'Keep-Alive')
        h.putheader('User-Agent', 'Bond/007 [en] (WinNT; U)')
        h.putheader('Content-Type', contentType)
        h.putheader('Content-Length', contentLength)
        h.endheaders()
        h.send(query)
        rcode, rmsg, headers= h.getreply()
        response = h.getfile()
        if return_type=='text':
            if rcode != 200:
                F='httpderror'
                msg = "error: %s, %s\n%s %s" % (rcode, self.uri, rmsg, response)
                tkMessageBox.showwarning("Could not contact pKD server",
                                         "Please check that your internet connection is functional.")
                return None
            else:
                return response.read()
        else:
            #
            # Return the file-like object
            #
            return response


