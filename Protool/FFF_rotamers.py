#
# Read the rotamer library for FFF
#
print 'Reading definition files for FFF'
import Protool, os, cPickle
dir_location=os.path.split(Protool.__file__)[0]
location=os.path.join(dir_location,'bbdep02.May.sortlib')
import FFF.FFFcontrol
Rotamerlib=FFF.FFFcontrol.Rotamer_class(location)
#
# Get the FFF aadef file
FFFaadef_file=os.path.split(FFF.FFFcontrol.__file__)[0]
FFFaadef_file=os.path.join(FFFaadef_file,'AA.DAT')

