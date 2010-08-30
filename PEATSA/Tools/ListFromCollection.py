#! /bin/env python
'''Recreates a mutation list file from a mutant collection'''
import PEAT_SA.Core as Core
import optparse, os, sys

#Parse args
usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-c", "--collection", dest="collection",
                  help="A PEAT-SA mutant collection directory.", metavar="COLLECTION")
		  
(options, args) = parser.parse_args()

if options.collection is None:
        print 'The location of a mutant collection must be provided'
        sys.exit(1)

components = os.path.split(options.collection)

#Do work
collection = Core.Data.MutantCollection(name=components[1], location=components[0])
mutations = collection.mutations()
list = Core.Data.MutationListFile(filename='mutationList', create=True)
for set in mutations:
	list.addMutant(set, autoUpdate=False)

list.writeToFile()
