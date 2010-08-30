#!/usr/bin/python -i

#
# Helper class for running FFF
#
import Protool, os
dir_location=os.path.split(Protool.__file__)[0]
location=os.path.join(dir_location,'bbdep02.May.sortlib')
text='Importing rotamer library for FFF from %s.....' %location
import FFFcontrol
Rotamerlib=FFFcontrol.Rotamer_class(location)
#
# Get the FFF aadef file
FFFaadef_file=os.path.split(FFFcontrol.__file__)[0]
FFFaadef_file=os.path.join(FFFaadef_file,'AA.DAT')


class FFF:

    def __init__(self,pdbfile):
        import FFFcontrol
        self.FFF=FFFcontrol.FFF()
        self.FFF.read_pdb(pdbfile)
        self.Model=FFFcontrol.model_class(self.FFF,Rotamerlib,FFFaadef_file)
        print 'Eenergy of res 133',self.Model.get_energy('E','133')
        print 'use X.Model to access model_class'
        return

if __name__=='__main__':
    X=FFF('impossible.pdb')
