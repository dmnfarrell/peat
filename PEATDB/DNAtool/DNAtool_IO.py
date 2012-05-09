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

import tkFileDialog, os
import tkMessageBox
import pickle
import PEATDB.DNA_sequence as DNA_sequence
import mutation

class DNA_IO:

    def project_open(self):
        """Open a project"""

        filename=tkFileDialog.askopenfilename(defaultextension='.DTP',
                                              initialdir=self.data['datadir'],
                                              parent=self.master,
                                              filetypes=[("DNATool project file","*.DTP"),
                                                         ("All files","*.*")])
        if filename:
            if os.path.isfile(filename):

                # Load the file
                try:
                    import pickle
                    fd=open(filename)
                    self.data=pickle.load(fd)
                    self.check_primers()
                    fd.close()
                    self.data['Project filename']=filename
                    self.open_pDB(overwrite=True)

                except:
                    import tkMessageBox
                    tkMessageBox.showwarning('Error reading file',
                                             'Please ensure that this is a valid DTP file',
                                             parent=self.master)
                    return

                self.assess_status()
                self.update_sequence_window()

            else:
                # File not found
                tkMessageBox.showwarning('File not found','Please select a valid file',
                                    parent=self.master)
        return

    def project_close(self):

        """Close the project"""
        if not self.data['Project saved']:
            if not tkMessageBox.askyesno("Project not saved", "Discard changes?",parent=self.master):
                return

        # Close everything
        self.master.destroy()
        self.__init__()
        return



    def project_saveas(self):
        # Get the filename

        filename=tkFileDialog.asksaveasfilename(defaultextension='.DTP',
                                                initialdir=self.data['datadir'],
                                                parent=self.master,
                                                filetypes=[("DNATool project file","*.DTP"),
                                                           ("All files","*.*")])
        if filename:
            self.data['Project filename']=filename
            self._save_project()
            self.assess_status()
        return

    def project_save(self):
        """Do we have a filename"""
        if self.data['Project filename']:
            self._save_project()
        else:
            tkMessageBox.showwarning('No project filename',
                                     'Please choose a project filename using Save As..',
                                     parent=self.master)
        return


    def _save_project(self):
        """Save everything including primer DB"""

        filename=self.data['Project filename']
        self.data['Project saved']=1
        print "Save data ", dir(self.data)
        #print self.data['primer_dict']
        fd=open(filename,'w')
        pickle.dump(self.data,fd)
        fd.close()

        self.assess_status()
        '''except:
            # We couldn't save??
            self.data['Project saved']=None
            import tkMessageBox
            tkMessageBox.showwarning('Save failed',
                                     'Could not save project. Make sure you have enough diskspace available'
                                     ,parent=self.master)'''

        return


    def open_file(self):
        """Open file dialog"""

        dnaseqfile=tkFileDialog.askopenfilename(defaultextension='.EAT',
                                                initialdir=self.data['datadir'],
                                                filetypes=[("PIR file","*.pir"),
                                                           ("FASTA file","*.txt"),
                                                           ("Clipped Seq","*.seq.clipped"),
                                                           ("GenBank file","*.gb"),
                                                           ("BSML file",".xml"),
                                                           ("All files","*.*")],
                                                parent=self.master)
        return dnaseqfile

    def dnaseq_read(self,newprotein=None,fileinput=None,variable='DNAseq',var2='dnaseqfile'):
        """Read in DNA Sequences and display them above the main sequence"""
        self.data['newprotein']=newprotein
        # If a file is provided use that, otherwise open the dialog
        if fileinput == None:
            dnaseqfile = self.open_file()
        else:
            dnaseqfile = fileinput
        self.data[var2]=dnaseqfile

        # Load the DNA sequence data
        #  - only simple file format supported at the moment

        if self.data[var2]:
            S=DNA_sequence.sequence()
            seqfile=self.data[var2]
            DNAseq=S.read_DNA_sequence(seqfile)

            # Check the DNA sequence
            ok,DNAseq=mutation.check_DNA(DNAseq)
            if not ok:
                # Give a warning
                import tkMessageBox
                tkMessageBox.showwarning("DNA sequence",
                                         'DNA sequence contains invalid characters: %s' %DNAseq
                                         ,parent=self.master)
                return
            self.data[variable]=DNAseq
            self.update_sequence_window()
            # Activate/Deactivate buttons
            self.assess_status()
        return dnaseqfile


    def save_dnaseq(self):
        """Save the DNA sequence"""
        print 'Not working yet'
        return


    def saveas_dnaseq(self):
        """Save the DNA sequence"""
        print 'Not working yet'
        return
