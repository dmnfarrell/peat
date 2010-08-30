#! /bin/env python
#
# Protein Engineering Analysis Tool Structure Analysis (PEATSA)
# Copyright (C) 2010 Michael Johnston & Jens Erik Nielsen
#
# Author: Michael Johnston
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
