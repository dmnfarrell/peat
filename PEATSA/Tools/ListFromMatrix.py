#! /usr/bin/python
'''Creates a mutation list from the Mutations column of a csv file

This script removes duplicates and provides error checking on the codes in the mutation column'''
import sys
import PEAT_SA.Core as Core

matrix = Core.Matrix.matrixFromCSVFile(sys.argv[1])
stream = Core.Data.MutationListFile(filename=sys.argv[2], mutationData="\n".join(matrix.mutations), create=True)
stream.removeDuplicates()
