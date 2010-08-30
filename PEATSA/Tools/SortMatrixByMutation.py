#! /usr/bin/env python
import PEATSA.Core as Core
import sys

def ShortMutationCodeCompare(code1, code2):

	#Codes are assumed to be extended
	data1 = Core.Utilities.ParseMutationCode(code1)
	data2 = Core.Utilities.ParseMutationCode(code2)

	if int(data1[1]) != int(data2[1]):
		return cmp(int(data1[1]), int(data2[1]))
	else:
		return cmp(data1[3], data2[3])


m = Core.Matrix.matrixFromCSVFile(sys.argv[1])

mutations = list(m.mutations)
mutations.sort(ShortMutationCodeCompare)

rows = []
for mutation in mutations:
	data = m.dataForMutation(mutation)
	rows.append(data)

sorted = Core.Matrix.PEATSAMatrix(rows=rows, headers=m.columnHeaders())

f = open('Sorted.csv', 'w+')
f.write(sorted.csvRepresentation())
f.close()
