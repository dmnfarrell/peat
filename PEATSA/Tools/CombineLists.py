#! /bin/env python
import sys
import PEAT_SA.Core as Core

files = sys.argv[1:]
print 'Combining ', files
mutationLists = [Core.Data.MutationListFile(file) for file in files]
newList = Core.Data.MutationListFile('newList', create=True)

for list in mutationLists:
	for set in list.mutantList():
		newList.addMutant(set, autoUpdate=False, ignoreDuplicates=True)

newList.removeDuplicates()
