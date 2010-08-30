'''Example script showing how to model a set of mutations on a pdb

Replace filename with an actual pdb'''
import sys, Protool
from PEAT_SA import Core

#Read in the pdb
filename = sys.argv[1]
pdb = Protool.structureIO()
pdb.readpdb(filename)

#Create the set of mutations to perform
mutationSet = Core.Data.MutationSet(name='MyClone1', code="A52H")

print '\n%s' % mutationSet
#To see the mutations relative to a sequence pass it a pdb
print 'Mutation I will perform is %s\n' % mutationSet.codeString(pdb)

#The MutantCollection class creates a set of mutants given a list of mutationSets
#You can 'load' any previously created set into another script aswell
#It can also be passed as an argument to ProteinDesignTool on the command line.
collection = Core.Data.MutantCollection(pdbFile=filename, name='ACollection', maxOverlap=1.4)
collection.createMutants(mutationList=[mutationSet]) 

print collection.mutantFiles()
