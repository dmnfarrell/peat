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

from Actions import DBActions

class ProteinEngineering:
    """This class handles general protein tasks and utilities, we should not
       put GUI stuff in here. Some items from old Administration.py are here"""
    
    @classmethod
    def alanineScan(self):
        """Perform a full Alanine scan on all the selected proteins
        Add all mutant proteins to the DB"""
        if not self.get_selected_proteins():
            return
        #
        # Open information window
        #
        P=GUI_helper.progress_window(self.master,'Alanine scan')
        P.update_progress(0.0)
        #
        # Loop over all proteins
        #
        DBi=self.data['DBinstance']
        DB=DBi.DB
        count=0
        for protein in self.sel_prots:
            #minres=DBi[protein]['ORF_selected']['aastart_number']
            #maxres=minres+DBi[protein]['ORF_selected']['length']
            sequence=DB[protein]['aaseq'][:]
            for resnum,restyp in sequence:
                print resnum,restyp
                new_name=DBi.suggest_name(protein,resnum,'ALA',restyp)
                mutation=resnum+':'+restyp+'='+'ALA'
                DBi.file_mutation_data(parent_name=protein,mutation=mutation,protein_name=new_name)
                #self.update_display()
                count=count+1
                P.update_progress(float(count)/float(len(self.sel_prots)*len(sequence)))
        P.close()
        #
        # Save the database
        #
        DBi.checkDB()
        DBi.saveDB()
        #
        # Update everything
        #
        self.update_all()
        return

    @classmethod
    def modelontheFly(self, DB, operations):
        """Produce the PDBlines for the mutant file"""
        
        self.MUT = DBActions.initProtool()

        # Create temporary array holding models
        if not getattr(DB,'tmp_3Dmodels',None):
            DB.tmp_3Dmodels={}
        # Get the structure from the parent
        parent = operations['parent']
        parentrec = DB.get(parent)
      
        if parentrec.hasStructure() == 'available':
            # Have we already modelled this one?
            if DB.tmp_3Dmodels.has_key(str(operations)):               
                pdblines = DB.tmp_3Dmodels[str(operations)]              
                import Protool
                X=Protool.structureIO()
                X.parsepdb(pdblines)
                return pdblines,X

            # Read parent structure
            pdblines = parentrec.Structure
            import Protool
            X=Protool.structureIO()
            X.parsepdb(pdblines)

            # Perform all the operations
            # First do all the deletes

            print operations['Rotamer_operations']
            for operation,uniqueid,atom in operations['Rotamer_operations']:
                #print operation,uniqueid,atom 
                if operation=='delete':                    
                    X.remove_atom(uniqueid)
            for operation,uniqueid,atom in operations['Rotamer_operations']:
                if operation=='add':                    
                    # Add atom                    
                    chainID=uniqueid.split(':')[0]
                    atomname=uniqueid.split(':')[-1]
                    residuenumber=uniqueid.split(':')[1]
                    resid='%s:%s' %(chainID,residuenumber)
                    status=X.add_atom(uniqueid=uniqueid,
                                      atomname=atomname,
                                      residuenumber=residuenumber,
                                      chainid=chainID,
                                      residuename=atom['RESNAME'],
                                      xcoord=atom['X'],ycoord=atom['Y'],zcoord=atom['Z'])
                    if not status:
                        self.parent.record_event('Modelling failed in model_on_the_fly for operations: %s' %(str(operations)))
                        return []
                    
                    # Update the arrays                    
                    X.Update()                    
                    # Rename the rest of the atoms                    
                    for this_atom in X.residues[resid]:
                        X.atoms[this_atom]['RESNAME']=atom['RESNAME']
                    X.Update()
                    #put mutation stuff into a tuple and return it too
                    #mutantinfo = (chainID,residuenumber)
            
            # Return the PDBlines and the instance            
            pdblines=X.writepdb('dummy',nowrite=1)
            
            # Store for next use            
            DB.tmp_3Dmodels[str(operations)]=pdblines[:]
            return pdblines, X
        else:
            print 'Parent does not have a structure'
            return [], None

