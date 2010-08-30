#! /usr/bin/env python
'''Converts a mutation list in the format ([ResidueNumber][MutationCode], ...) e.g. 54A,63W,84A

into a PEAT_SA format - assuming a homodimer.

Usage:
	ConvertHIVDBList.py [inputListFile] [outputListFile]
'''
import sys
import PEAT_SA.Core as Core

f = open(sys.argv[1])
originalLines = f.readlines()
originalLines = [line.strip('\n') for line in originalLines]
f.close()

lines = [line.split(',') for line in originalLines]
listStream = Core.Data.MutationListFile(sys.argv[2], create=True)
sets = []

print 'Parsing HIVDB list'
for line in lines:
	chainA = ['A'+el for el in line]
	chainB = ['B'+el for el in line]
	chainA = "+".join(chainA)
	chainB = "+".join(chainB)
	code = chainA+'+'+chainB
	mutationSet = Core.Data.MutationSet(code)
	sets.append(mutationSet)

#Create matrix of peatsa code HIVDB code pairs
codes = ["+".join(mutationSet.reducedMutationCodes()) for mutationSet in sets]
lengthCol = [len(mutationSet.mutationCodes()) for mutationSet in sets]
#Remove commas since they will be seen as entry delimiters in a csv file
originalLines = [line.replace(',', '+') for line in originalLines]
rows = zip(codes, originalLines)
matrix = Core.Matrix.Matrix(rows=rows, headers=['PEATSA', 'HIVDB'])
matrix.addColumn(lengthCol)
matrix.setColumnHeader(2, 'Length')

print 'Writing HIV->PEAT_SA map to Map.csv'
stream = open('Map.csv', 'w')
stream.write(matrix.csvRepresentation())
stream.close()

#Create mutation list
sets.sort(lambda x,y: cmp(len(x.mutationCodes()), len(y.mutationCodes())))
for mutationSet in sets:
	listStream.addMutant(mutationSet, autoUpdate=False)

print 'Writing peatsa mutation list to %s' % sys.argv[2]
listStream.writeToFile()

#All codes
codes = []
[codes.extend(line) for line in lines]
nonRedundant = list(set(codes))
for code in nonRedundant:
	print code






